#include "sdf_scanner/sdf_scan.hpp"
#include "duckdb/common/allocator.hpp"
#include "duckdb/common/helper.hpp"
#include "duckdb/common/multi_file_reader.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/unique_ptr.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
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

    vector<string> cur_row;
    auto cur_mol = gstate.mol_supplier->next();
    //! Go through each column specified and store the property in a vector
    //! This represents one row in the "table".
    for (idx_t i = 0; i < bind_data.names.size(); i++) {
      //! NOTE: using the RDKit MolSupplier, if the record cannot be parsed, it
      //! is because the molecule cannot be parsed. The other columns are
      //! probably not null, but right now nothing of that record is
      //! stored because that cannot be accessed with the way the MolSupplier
      //! is implemented.
      //! TODO: write a parser to parse the non-molecule parts of the SDF.
      //! This may enable filter pushdown and also allow the other fields
      //! to be returned even if the molecule cannot be parsed
      if (cur_mol) {
        std::string prop;
        cur_mol->getPropIfPresent(bind_data.names[i], prop);
        cur_row.emplace_back(prop);
      } else {
        cur_row.emplace_back("");
      }
    }
    lstate.rows.emplace_back(cur_row);
    lstate.scan_count++;
    gstate.offset++;
  }
}

SDFLocalTableFunctionState::SDFLocalTableFunctionState(
    ClientContext &context_p, SDFScanGlobalState &gstate_p)
    : state(context_p, gstate_p) {}

//
// unique_ptr<LocalTableFunctionState>
// SDFLocalTableFunctionState::Init(ExecutionContext &context_p,
//                                  TableFunctionInitInput &input_p,
//                                  GlobalTableFunctionState *global_state_p)
//                                  {
//   auto &gstate = global_state_p->Cast<SDFGlobalTableFunctionState>();
//   auto result =
//       make_uniq<SDFLocalTableFunctionState>(context_p.client,
//       gstate.state);
//   return std::move(result);
// }
//
// idx_t SDFLocalTableFunctionState::GetBatchIndex() const {
//   return state.batch_index;
// }
//
} // namespace duckdb
