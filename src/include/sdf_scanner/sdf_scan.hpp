#pragma once
#include "GraphMol/FileParsers/MolSupplier.h"
#include "duckdb/common/multi_file_list.hpp"
#include "duckdb/common/multi_file_reader.hpp"
#include "duckdb/common/multi_file_reader_options.hpp"
#include "duckdb/common/unique_ptr.hpp"
#include "duckdb/execution/execution_context.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/function_set.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "sdf_scanner/buffered_sdf_reader.hpp"
#include "sdf_scanner/plain_sdf_reader.hpp"
namespace duckdb {

//! The code here is based off the duckdb JSON extension

struct SDFScanData : public TableFunctionData {
public:
  SDFScanData();

  void Bind(ClientContext &context, TableFunctionBindInput &input);

  //! Multi-file reader stuff
  MultiFileReaderBindData reader_bind;
  //! The files we're reading
  vector<string> files;
  //! Multi-file reader options
  MultiFileReaderOptions file_options;
  //! All column names (in order) specified by the query for projection
  vector<string> names;
  idx_t maximum_object_size = 16777216;
};

struct SDFScanGlobalState {
public:
  SDFScanGlobalState(ClientContext &context, const SDFScanData &bind_data);

public:
  //! Bound data
  const SDFScanData &bind_data;
  unique_ptr<RDKit::v2::FileParsers::SDMolSupplier> mol_supplier;

  //! Column names that we're actually reading (after projection pushdown)
  vector<string> names;
};

struct SDFScanLocalState {
public:
  SDFScanLocalState(ClientContext &context, SDFScanGlobalState &gstate);

public:
  //! Retrieves the next chunk of SDF records from the global state
  //! and populates the local state with its properties.
  //! If a molecule is available, its properties are extracted and stored in
  //! `lstate.properties`, and `lstate.scan_count` is incremented. If no
  //! molecule is available, logs a failure message and leaves
  //! `lstate.scan_count` at zero.
  void ExtractNextChunk(SDFScanGlobalState &gstate, SDFScanLocalState &lstate,
                        SDFScanData &bind_data);

public:
  //! The number of records successfully scanned from the SDF.
  //! This is used to indicate duckdb the number of rows
  //! returned by the scanning function. If this value is zero, duckdb will then
  //! know to not call the function again. That will end the scan.
  idx_t scan_count;

  //! The properties for each record in the current chunk of records being
  //! scanned. Each record is a "row", each property is a "column"
  vector<vector<string>> rows;

private:
  //! Bind data
  const SDFScanData &bind_data;
};

struct SDFGlobalTableFunctionState : public GlobalTableFunctionState {
public:
  SDFGlobalTableFunctionState(ClientContext &context,
                              TableFunctionInitInput &input);
  static unique_ptr<GlobalTableFunctionState>
  Init(ClientContext &context, TableFunctionInitInput &input);

public:
  SDFScanGlobalState state;
};

struct SDFLocalTableFunctionState : public LocalTableFunctionState {
public:
  SDFLocalTableFunctionState(ClientContext &context,
                             SDFScanGlobalState &gstate);
  static unique_ptr<LocalTableFunctionState>
  Init(ExecutionContext &context, TableFunctionInitInput &input,
       GlobalTableFunctionState *global_state);
  idx_t GetBatchIndex() const;

public:
  SDFScanLocalState state;
};

class SDFScanFunction {
public:
  static TableFunctionSet GetFunctionSet();

  static void SDFScanImplementation(ClientContext &context,
                                    TableFunctionInput &data_p,
                                    DataChunk &output);

  static unique_ptr<FunctionData> SDFScanBind(ClientContext &context,
                                              TableFunctionBindInput &input,
                                              vector<LogicalType> &return_types,
                                              vector<string> &names);

  static unique_ptr<GlobalTableFunctionState>
  SDFScanInitGlobal(ClientContext &context, TableFunctionInitInput &input);

  static unique_ptr<LocalTableFunctionState>
  SDFScanInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                   GlobalTableFunctionState *gstate);

private:
};

struct SDFReadBindData : public TableFunctionData {
  shared_ptr<MultiFileList> file_list;
  unique_ptr<MultiFileReader> multi_file_reader;
  vector<string> names;
  vector<LogicalType> types;
};

struct SDFReadGlobalState : public GlobalTableFunctionState {
  explicit SDFReadGlobalState();
  mutex lock;
  vector<idx_t> projection_ids;
  optional_ptr<TableFilterSet> filters;
  //! The file list to scan
  MultiFileList &file_list;
  //! The scan over the file_list
  MultiFileListScanData file_list_scan;
  unique_ptr<MultiFileReaderGlobalState> multi_file_reader_state;
};

} // namespace duckdb
