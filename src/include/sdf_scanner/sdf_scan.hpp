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

  //! Buffer manager allocator
  // Allocator &allocator;
  //! The current buffer capacity
  // idx_t buffer_capacity;

  //! Column names that we're actually reading (after projection pushdown)
  vector<string> names;
  vector<column_t> column_ids;
  // vector<ColumnIndex> column_indices;

  mutex lock;
  //! One SDF reader per file
  // vector<optional_ptr<PlainSDFReader>> sdf_readers;
  //! Current file/batch index
  // atomic<idx_t> file_index;
  // atomic<idx_t> batch_index;

  //! Current number of threads active
  // idx_t system_threads;
  //! Whether we enable parallel scans (only if less files than threads)
  // bool enable_parallel_scans;
};

struct SDFScanLocalState {
public:
  SDFScanLocalState(ClientContext &context, SDFScanGlobalState &gstate);

public:
  void ReadNext(SDFScanGlobalState &gstate);

  const MultiFileReaderData &GetReaderData() const;

public:
  //! Current scan data
  idx_t scan_count;
  //! The current molecule that has been read by the mol_supplier
  unique_ptr<RDKit::ROMol> cur_mol;

  //! The properties of the current record
  vector<string> properties;
  // string units[STANDARD_VECTOR_SIZE];

  //! Batch index for order-preserving parallelism
  // idx_t batch_index;

  //! For determining average tuple size
  // idx_t total_read_size;
  // idx_t total_tuple_count;

private:
  AllocatedData AllocateBuffer(SDFScanGlobalState &gstate);

  void ParseNextChunk(SDFScanGlobalState &gstate);

  void ParseSDF(const idx_t sdf_start, const idx_t sdf_size,
                const idx_t remaining);
  void ThrowObjectSizeError(const idx_t object_size);

  //! Must hold the lock
  void TryIncrementFileIndex(SDFScanGlobalState &gstate) const;
  bool IsParallel(SDFScanGlobalState &gstate) const;

private:
  //! Bind data
  const SDFScanData &bind_data;
  //! Thread-local allocator
  // Allocator &allocator;

  //! Current reader and buffer handle
  // optional_ptr<PlainSDFReader> current_reader;
  //! Whether this is the last batch of the file
  // bool is_last;

  //! The current main filesystem
  // FileSystem &fs;

  //! For some filesystems (e.g. S3), using a filehandle per thread increases
  //! performance
  // unique_ptr<FileHandle> thread_local_filehandle;

  //! Current buffer read info
  // char *buffer_ptr;
  // idx_t buffer_size;
  // idx_t buffer_offset;
  // idx_t prev_buffer_remainder;
  // idx_t lines_or_objects_in_buffer;

  //! Buffer to reconstruct split values
  // AllocatedData reconstruct_buffer;
};

struct SDFGlobalTableFunctionState : public GlobalTableFunctionState {
public:
  SDFGlobalTableFunctionState(ClientContext &context,
                              TableFunctionInitInput &input);
  static unique_ptr<GlobalTableFunctionState>
  Init(ClientContext &context, TableFunctionInitInput &input);
  // idx_t MaxThreads() const override;

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
