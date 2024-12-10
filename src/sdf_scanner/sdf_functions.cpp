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

void ExtractNextChunk(SDFScanGlobalState &gstate, SDFScanLocalState &lstate,
                      SDFScanData &bind_data) {
  auto cur = gstate.mol_supplier->next();
  //! this holds the number of records scanned
  //! reset to zero each time this function is called
  //! If nothing gets scanned, the scan_count will be just zero
  lstate.scan_count = 0;
  lstate.properties.clear();
  if (cur) {
    for (idx_t i = 0; i < bind_data.names.size(); i++) {
      std::string prop;
      std::cout << bind_data.names[i] << std::endl;
      cur->getPropIfPresent(bind_data.names[i], prop);
      std::cout << prop << std::endl;
      lstate.properties.emplace_back(prop);
      std::cout << lstate.properties[i] << std::endl;
    }
    lstate.scan_count++;
  } else {
    std::cout << "could not convert "
              << gstate.mol_supplier->getItemText(lstate.scan_count)
              << std::endl;
  }
}

static void ReadSDFFunction(ClientContext &context, TableFunctionInput &data_p,
                            DataChunk &output) {

  auto &gstate = data_p.global_state->Cast<SDFGlobalTableFunctionState>().state;
  auto &lstate = data_p.local_state->Cast<SDFLocalTableFunctionState>().state;
  auto bind_data = data_p.bind_data->Cast<SDFScanData>();
  // auto offset = lstate.scan_count;
  // auto count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, 10 - offset);
  ExtractNextChunk(gstate, lstate, bind_data);
  std::cout << "COUNT: " << lstate.scan_count << std::endl;

  // set to the number of rows returned
  output.SetCardinality(lstate.scan_count);

  // If SetCardinality is > 0, it goes in an infinite loop
  // mol_supplier.close();
  std::cout << "number of tuples: " << output.size() << std::endl;
  std::cout << "number of columns: " << output.data.size() << std::endl;
  if (lstate.properties.size() > 0) {

    for (idx_t i = 0; i < bind_data.names.size(); i++) {
      auto k = lstate.properties[i];
      std::cout << "bind data " << i << ": " << bind_data.names[i] << std::endl;
      std::cout << "PROPERTY " << i << ": " << k << std::endl;
      output.SetValue(i, 0, Value(k));
    }
  }
  std::cout << output.ToString() << std::endl;
  std::cout << output.ColumnCount() << std::endl;
  /// NOTE: everything is just string_t right now
  /// make sure to declare only VARCHAR in COLUMNS
  /// Perhaps, don't even make it an option for the first iteration
  /// Just make it default all VARCHAR. So user could just specify
  /// COLUMNS=[names, of, columns] ?
  ///
  // auto data = FlatVector::GetData<string_t>(output.data[0]);
  // auto &validity = FlatVector::Validity(output.data[0]);
  //
  // std::cout << data->GetString() << std::endl;
  // FlatVector::SetData(output.data[1], output.data[1].GetData());
  // auto data1 = FlatVector::GetData<string_t>(output.data[1]);
  // auto &validity1 = FlatVector::Validity(output.data[1]);
  //
  // std::cout << data1->GetString() << std::endl;

  // Left off at: why it doesn't return as rows in the query?

  // MultiFileConstantEntry entry(0, output.data[0].GetValue(0));
  // MultiFileConstantEntry entry1(1, output.data[1].GetValue(0));
  // vector<MultiFileConstantEntry> constant_map;
  // constant_map.emplace_back(entry);
  // constant_map.emplace_back(entry1);
  //
  // MultiFileReaderData reader_data;
  // reader_data.constant_map = constant_map;
  // MultiFileReaderBindData reader_bind;
  // if (output.size() != 0) {
  //   MultiFileReader().FinalizeChunk(context, reader_bind, reader_data,
  //   output,
  //                                   nullptr);
  // }
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
  //! TODO: should I use the multi file reader or use the SDF Mol Supplier?
  // MultiFileReader().BindOptions(bind_data->file_options, file_list,
  //                               return_types, names,
  //                               bind_data->reader_bind);
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
