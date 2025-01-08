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
namespace duckdb {

//! The code here is based off the duckdb JSON extension

struct SDFScanData : public TableFunctionData {
public:
  SDFScanData();

  void Bind(ClientContext &context, TableFunctionBindInput &input);

  //! The files we're reading
  vector<string> files;
  //! All column names (in order) specified by the query for projection
  vector<string> names;
  //! All column types (in order) specified by the query for projection
  //!
  //! This is not combined with the names field because the TableFunction
  //! constructor expects a vector<string> for the column names field
  //! specifically
  vector<string> types;

  //! If a Mol type is requested, it is necessary to know what index it is
  //! in the column vector in order to handle the data differently during
  //! read of the SDF
  //! A value of -1 means there is no mol_col_idx set because Mol type was
  //! not requested
  short mol_col_idx = -1;
};

struct SDFScanGlobalState {
public:
  SDFScanGlobalState(ClientContext &context, const SDFScanData &bind_data);

public:
  //! Bound data
  const SDFScanData &bind_data;
  unique_ptr<RDKit::v2::FileParsers::SDMolSupplier> mol_supplier;
  //! the number of records in the SDF file
  idx_t length;
  //! The record number of the last record we have scanned in the SDF file;
  //! where we are in the SDF
  idx_t offset;

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
  //! the local state. The number of records scanned is also incremented
  //! accordingly, and stored in the local state.
  void ExtractNextChunk(SDFScanGlobalState &gstate, SDFScanLocalState &lstate,
                        SDFScanData &bind_data);

public:
  //! The number of records successfully scanned from the SDF.
  //! This is used to indicate to duckdb the number of rows
  //! returned by the scanning function. If this value is zero, duckdb will then
  //! know to not call the function again. That will end the scan.
  idx_t scan_count;

  //! The properties for each record in the current chunk of records being
  //! scanned. Each record is a "row", each property is a "column". The outer
  //! vector is a vector of "rows" (the inner vector). Each entry in the inner
  //! vector is a "column", or property extracted from the SDF fields
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

public:
  SDFScanLocalState state;
};

struct SDFScan {
public:
  static double ScanProgress(ClientContext &context,
                             const FunctionData *bind_data_p,
                             const GlobalTableFunctionState *global_state);
};

} // namespace duckdb
