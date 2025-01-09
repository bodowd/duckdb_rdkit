#include "sdf_scanner/sdf_functions.hpp"
#include "duckdb/common/assert.hpp"
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
#include "types.hpp"
#include <iterator>

namespace duckdb {

static void ReadSDFFunction(ClientContext &context, TableFunctionInput &data_p,
                            DataChunk &output) {

  auto &gstate = data_p.global_state->Cast<SDFGlobalTableFunctionState>().state;
  auto &lstate = data_p.local_state->Cast<SDFLocalTableFunctionState>().state;
  auto bind_data = data_p.bind_data->Cast<SDFScanData>();

  lstate.ExtractNextChunk(gstate, lstate, bind_data);

  D_ASSERT(lstate.scan_count == lstate.rows.size());
  //! set to the number of rows scanned
  //! If the cardinality is zero, it will signal to duckdb to not run the read
  //! function anymore because the scan is done. This is not equal to the number
  //! of rows returned, if perhaps during predicate pushdown, no rows match.
  //! If we are not done scanning the file, we would need to keep scanning
  //! to see if there are other records that match the predicate
  output.SetCardinality(lstate.scan_count);
  if (bind_data.mol_col_idx > -1) {
    //! The mol_col is a reference to the DataChunk vector
    //! This vector is then converted to a FlatVector with the string_t type
    //! so that we can add blobs to it, which may have invalid UTF8
    auto &mol_col = output.data[bind_data.mol_col_idx];
    auto col_data = FlatVector::GetData<string_t>(mol_col);
    for (idx_t i = 0; i < lstate.rows.size(); i++) {
      for (idx_t j = 0; j < bind_data.names.size(); j++) {
        if (bind_data.mol_col_idx > -1 && j == bind_data.mol_col_idx) {
          auto val = string_t(lstate.rows[i][j]);
          col_data[i] = StringVector::AddStringOrBlob(mol_col, val);
        }
      }
    }
  }

  //! For each record/row scanned, set the value of each
  //! column in the output DataChunk
  for (idx_t i = 0; i < lstate.rows.size(); i++) {
    for (idx_t j = 0; j < bind_data.names.size(); j++) {
      auto val = lstate.rows[i][j];
      if (bind_data.mol_col_idx > -1 && j != bind_data.mol_col_idx) {
        if (val == "") {
          output.SetValue(j, i, Value(nullptr));
        } else {
          output.SetValue(j, i, Value(val));
        }
      }
    }
  }
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
      vector<string> types;
      D_ASSERT(StructType::GetChildCount(child_type) == struct_children.size());
      for (idx_t i = 0; i < struct_children.size(); i++) {
        auto &name = StructType::GetChildName(child_type, i);
        auto &type = struct_children[i];
        names.push_back(name);
        types.push_back(type.ToString());
        if (type.type().id() != LogicalTypeId::VARCHAR) {
          throw BinderException("read_sdf \"columns\" parameter type "
                                "specification must be VARCHAR.");
        }

        if (type.ToString() == duckdb_rdkit::Mol().ToString()) {
          bind_data->mol_col_idx = i;
          return_types.emplace_back(duckdb_rdkit::Mol());
        } else {
          //! All columns that are not Mol type should be converted to a normal
          //! duckdb LogicalType
          return_types.emplace_back(
              TransformStringToLogicalType(StringValue::Get(type), context));
        }
      }
      D_ASSERT(names.size() == return_types.size());
      if (names.empty()) {
        throw BinderException(
            "read_sdf \"columns\" parameter needs at least one column.");
      }
      bind_data->names = names;
      bind_data->types = types;
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

  functions.push_back(GetReadSDFTableFunction());

  return functions;
}

TableFunctionSet SDFFunctions::GetReadSDFTableFunction() {
  TableFunction table_function({LogicalType::VARCHAR}, ReadSDFFunction,
                               ReadSDFBind, SDFGlobalTableFunctionState::Init,
                               SDFLocalTableFunctionState::Init);
  table_function.name = "read_sdf";
  table_function.named_parameters["columns"] = LogicalType::ANY;
  table_function.table_scan_progress = SDFScan::ScanProgress;
  return MultiFileReader::CreateFunctionSet(table_function);
}

} // namespace duckdb
