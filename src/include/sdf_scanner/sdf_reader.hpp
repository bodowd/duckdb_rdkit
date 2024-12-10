// #pragma once
//
// #include "duckdb/common/file_system.hpp"
// #include "duckdb/common/types/data_chunk.hpp"
// #include "duckdb/common/vector.hpp"
// #include "duckdb/function/function_set.hpp"
// #include "duckdb/main/client_context.hpp"
// namespace duckdb {
//
// struct SDFReaderScanState {
//   unique_ptr<FileHandle> file_handle;
//   bool finished;
// };
//
// class SDFReader {
// public:
//   SDFReader(ClientContext &context, string filename);
//   ~SDFReader();
//
//   void InitializeScan(ClientContext &context, SDFReaderScanState &state);
//   void Scan(SDFReaderScanState &state, DataChunk &result);
//   bool ScanInternal(SDFReaderScanState &state, DataChunk &output);
//
//   FileSystem &fs;
//   string file_name;
//
// private:
//   unique_ptr<FileHandle> file_handle;
// };
// } // namespace duckdb
