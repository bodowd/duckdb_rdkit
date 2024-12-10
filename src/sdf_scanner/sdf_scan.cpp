#include "sdf_scanner/sdf_scan.hpp"
#include "duckdb/common/allocator.hpp"
#include "duckdb/common/helper.hpp"
#include "duckdb/common/multi_file_reader.hpp"
#include "duckdb/common/unique_ptr.hpp"
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
}

SDFScanLocalState::SDFScanLocalState(ClientContext &context_p,
                                     SDFScanGlobalState &gstate_p)
    : scan_count(0), bind_data(gstate_p.bind_data) {}

SDFGlobalTableFunctionState::SDFGlobalTableFunctionState(
    ClientContext &context, TableFunctionInitInput &input)
    : state(context, input.bind_data->Cast<SDFScanData>()) {}

void SDFScanLocalState::ReadNext(SDFScanGlobalState &gstate) {
  cur_mol = gstate.mol_supplier->next();
  std::string tmp;
  cur_mol->getPropIfPresent("ChEBI Name", tmp);

  std::cout << "!!!!!!!! In read next: " << tmp << std::endl;
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
