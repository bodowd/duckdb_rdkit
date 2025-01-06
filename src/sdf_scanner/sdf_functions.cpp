#include "sdf_scanner/sdf_functions.hpp"
#include "duckdb/common/enums/vector_type.hpp"
#include "duckdb/common/multi_file_list.hpp"
#include "duckdb/common/multi_file_reader.hpp"
#include "duckdb/common/operator/cast_operators.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "mol_formats.hpp"
#include "sdf_scanner/sdf_scan.hpp"

namespace duckdb {

static void ReadSDFFunction(ClientContext &context, TableFunctionInput &data_p,
                            DataChunk &output) {

  auto &gstate = data_p.global_state->Cast<SDFGlobalTableFunctionState>().state;
  auto &lstate = data_p.local_state->Cast<SDFLocalTableFunctionState>().state;
  auto bind_data = data_p.bind_data->Cast<SDFScanData>();

  lstate.ExtractNextChunk(gstate, lstate, bind_data);

  D_ASSERT(lstate.scan_count == lstate.rows.size());
  // set to the number of rows returned
  output.SetCardinality(lstate.scan_count);

  //! For each record scanned, set the value of each
  //! column in the output DataChunk
  for (idx_t i = 0; i < lstate.rows.size(); i++) {
    for (idx_t j = 0; j < bind_data.names.size(); j++) {
      auto val = lstate.rows[i][j];
      if (val == "") {
        output.SetValue(j, i, Value(nullptr));
      } else {
        output.SetValue(j, i, Value(val));
      }
    }
  }
  /// NOTE: everything is just string_t right now
  /// make sure to declare only VARCHAR in COLUMNS
  /// Perhaps, don't even make it an option for the first iteration
  /// Just make it default all VARCHAR. So user could just specify
  /// COLUMNS=[names, of, columns] ?
}

unique_ptr<FunctionData> ReadSDFBind(ClientContext &context,
                                     TableFunctionBindInput &input,
                                     vector<LogicalType> &return_types,
                                     vector<string> &names) {
  auto bind_data = make_uniq<SDFScanData>();
  bind_data->Bind(context, input);
  for (auto &kv : input.named_parameters) {
    //! kv.first represents the parameter name
    if (kv.second.IsNull()) {
      throw BinderException("Cannot use NULL as function argument");
    }
    auto loption = StringUtil::Lower(kv.first);
    if (kv.second.IsNull()) {
      throw BinderException("read_sdf parameter \"%s\" cannot be NULL.",
                            loption);
    }
    if (loption == "columns") {
      //! kv.second is the struct of columns
      //! Its shape is {name_of_column: logical_type_of_column}
      auto &child_type = kv.second.type();
      if (child_type.id() != LogicalTypeId::STRUCT) {
        throw BinderException(
            "read_sdf \"columns\" parameter requires a struct as input.");
      }
      //! struct_children are the types of the columns
      auto &struct_children = StructValue::GetChildren(kv.second);
      D_ASSERT(StructType::GetChildCount(child_type) == struct_children.size());
      for (idx_t i = 0; i < struct_children.size(); i++) {
        auto &name = StructType::GetChildName(child_type, i);
        auto &val = struct_children[i];
        names.push_back(name);
        if (val.type().id() != LogicalTypeId::VARCHAR) {
          throw BinderException("read_sdf \"columns\" parameter type "
                                "specification must be VARCHAR.");
        }
        return_types.emplace_back(
            TransformStringToLogicalType(StringValue::Get(val), context));
      }
      D_ASSERT(names.size() == return_types.size());
      if (names.empty()) {
        throw BinderException(
            "read_sdf \"columns\" parameter needs at least one column.");
      }
      bind_data->names = names;
    }

    //! TODO: not implemented yet
    // else if (loption == "convert_strings_to_integers") {
    // bind_data->convert_strings_to_integers =
    // BooleanValue::Get(kv.second);
    // }
  }

  //! get the files
  SimpleMultiFileList file_list(std::move(bind_data->files));
  bind_data->files = file_list.GetAllFiles();
  if (bind_data->files.size() > 1) {
    throw NotImplementedException(
        "Reading more than one sdf file is currently not supported.");
  }

  return std::move(bind_data);
}

vector<TableFunctionSet> SDFFunctions::GetTableFunctions() {
  vector<TableFunctionSet> functions;

  functions.push_back(GetReadSDFFunction());

  return functions;
}

TableFunction GetReadSDFTableFunction() {
  TableFunction table_function({LogicalType::VARCHAR}, ReadSDFFunction,
                               ReadSDFBind, SDFGlobalTableFunctionState::Init,
                               SDFLocalTableFunctionState::Init);
  table_function.name = "read_sdf";
  table_function.named_parameters["columns"] = LogicalType::ANY;

  return table_function;
}

TableFunctionSet SDFFunctions::GetReadSDFFunction() {

  auto table_function = GetReadSDFTableFunction();
  return MultiFileReader::CreateFunctionSet(table_function);
}

} // namespace duckdb
