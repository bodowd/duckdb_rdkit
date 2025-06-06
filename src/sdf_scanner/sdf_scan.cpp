#include "sdf_scanner/sdf_scan.hpp"
#include "duckdb/common/allocator.hpp"
#include "duckdb/common/helper.hpp"
#include "duckdb/common/multi_file_reader.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/unique_ptr.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "mol_formats.hpp"
#include "types.hpp"
#include "umbra_mol.hpp"
#include <memory>

namespace duckdb {
SDFScanData::SDFScanData() {}

void SDFScanData::Bind(ClientContext &context, TableFunctionBindInput &input) {
  auto multi_file_reader = MultiFileReader::Create(input.table_function);
  auto file_list = multi_file_reader->CreateFileList(context, input.inputs[0]);

  files = file_list->GetAllFiles();
}

unique_ptr<LocalTableFunctionState>
SDFLocalTableFunctionState::Init(ExecutionContext &context,
                                 TableFunctionInitInput &,
                                 GlobalTableFunctionState *global_state) {
  auto &gstate = global_state->Cast<SDFGlobalTableFunctionState>();
  auto result =
      make_uniq<SDFLocalTableFunctionState>(context.client, gstate.state);
  return std::move(result);
}

unique_ptr<GlobalTableFunctionState>
SDFGlobalTableFunctionState::Init(ClientContext &context,
                                  TableFunctionInitInput &input) {
  auto &bind_data = input.bind_data->Cast<SDFScanData>();
  auto result = make_uniq<SDFGlobalTableFunctionState>(context, input);
  auto &gstate = result->state;

  return std::move(result);
}

SDFScanGlobalState::SDFScanGlobalState(ClientContext &context_p,
                                       const SDFScanData &bind_data_p)
    : bind_data(bind_data_p) {
  //! open the sd file
  mol_supplier =
      make_uniq<RDKit::v2::FileParsers::SDMolSupplier>(bind_data.files[0]);
  length = mol_supplier->length();
  offset = 0;
}

SDFScanLocalState::SDFScanLocalState(ClientContext &context_p,
                                     SDFScanGlobalState &gstate_p)
    : scan_count(0), bind_data(gstate_p.bind_data) {}

SDFGlobalTableFunctionState::SDFGlobalTableFunctionState(
    ClientContext &context, TableFunctionInitInput &input)
    : state(context, input.bind_data->Cast<SDFScanData>()) {}

void SDFScanLocalState::ExtractNextChunk(SDFScanGlobalState &gstate,
                                         SDFScanLocalState &lstate,
                                         SDFScanData &bind_data) {

  //! This holds the number of records scanned in this current
  //! function call.
  //! Reset to zero each time this function is called.
  //! If nothing gets scanned, the scan_count will be just zero
  //! and duckdb will be signalled that the scanning is complete
  lstate.scan_count = 0;
  lstate.rows.clear();

  while (lstate.scan_count < STANDARD_VECTOR_SIZE &&
         !gstate.mol_supplier->atEnd()) {
    auto remaining = gstate.length - lstate.scan_count;
    if (remaining == 0) {
      break;
    }

    bool printed_warning = false;
    vector<string> cur_row;
    auto cur_mol = gstate.mol_supplier->next();
    //! Go through each column specified and store the property in a vector
    //! This represents one row in the "table".
    for (idx_t i = 0; i < bind_data.names.size(); i++) {
      //! NOTE: using the RDKit MolSupplier, if the record cannot be parsed,
      //! it is because the molecule cannot be parsed. The other columns are
      //! probably not null, but right now nothing of that record is
      //! stored because that cannot be accessed with the way the MolSupplier
      //! is implemented.
      //!
      //! TODO: write a parser to parse the non-molecule parts of the SDF.
      //! This may enable filter pushdown and also allow the other fields
      //! to be returned even if the molecule cannot be parsed
      if (cur_mol) {
        //! The column at the current iteration of i is the Mol type
        //! In this case, we should convert the molecule object
        //! to the "umbra" mol in duckdb_rdkit
        if (bind_data.types[i] == duckdb_rdkit::Mol().ToString()) {
          auto res = duckdb_rdkit::get_umbra_mol_string(*cur_mol);
          cur_row.emplace_back(res);
        } else {
          //! Otherwise, it is a normal property column
          std::string prop;
          cur_mol->getPropIfPresent(bind_data.names[i], prop);
          cur_row.emplace_back(prop);
        }
      } else {
        //! only print the warning once per record, not once per column
        //! of the record
        if (!printed_warning) {
          std::cout << "Molecule could not be constructed at record #: "
                    << gstate.offset << std::endl;
          printed_warning = true;
        }
        cur_row.emplace_back("");
      }
    }
    lstate.rows.emplace_back(cur_row);
    lstate.scan_count++;
    gstate.offset++;
  }
}

void SDFScan::AutoDetect(ClientContext &context, SDFScanData &bind_data,
                         vector<LogicalType> &return_types,
                         vector<string> &names) {
  //! open the sd file to scan the first record
  auto mol_supplier =
      make_uniq<RDKit::v2::FileParsers::SDMolSupplier>(bind_data.files[0]);
  while (!mol_supplier->atEnd()) {

    auto cur_mol = mol_supplier->next();
    std::set<string> seen;
    if (cur_mol) {
      for (auto p : cur_mol->getPropList()) {
        //! These are props seem to be added by RDKit...these are there
        //! even if the SDF doesn't contain these properties
        if (p != "__computedProps" && p != "_Name" && p != "_MolFileInfo" &&
            p != "_MolFileComments" && p != "_MolFileChiralFlag" &&
            p != "numArom" && p != "_StereochemDone") {
          names.push_back(p);
          bind_data.types.emplace_back(
              LogicalTypeIdToString(LogicalTypeId::VARCHAR));
          return_types.emplace_back(TransformStringToLogicalType(
              StringValue::Get("VARCHAR"), context));
        }
      }
      break;
    }
  }
  mol_supplier->close();

  //! The molecule block is not in the getPropList
  names.push_back("mol");
  bind_data.types.push_back("Mol");
  bind_data.mol_col_idx = names.size() - 1;
  return_types.emplace_back(duckdb_rdkit::Mol());
  bind_data.names = names;
}

SDFLocalTableFunctionState::SDFLocalTableFunctionState(
    ClientContext &context_p, SDFScanGlobalState &gstate_p)
    : state(context_p, gstate_p) {}

double SDFScan::ScanProgress(ClientContext &, const FunctionData *,
                             const GlobalTableFunctionState *global_state) {
  auto &gstate = global_state->Cast<SDFGlobalTableFunctionState>().state;
  return 100.0 * (double)gstate.offset / (double)gstate.length;
}

} // namespace duckdb
