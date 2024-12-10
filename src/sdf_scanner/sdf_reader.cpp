// #include "sdf_scanner/sdf_reader.hpp"
// #include "duckdb/common/file_open_flags.hpp"
// #include "duckdb/common/file_system.hpp"
// #include "duckdb/main/client_context.hpp"
// #include <iostream>
// namespace duckdb {
//
// SDFReader::SDFReader(ClientContext &context, string filename)
//     : fs(FileSystem::GetFileSystem(context)) {
//   file_name = std::move(filename);
//   file_handle = fs.OpenFile(file_name, FileFlags::FILE_FLAGS_READ);
//   if (!file_handle->CanSeek()) {
//     throw NotImplementedException(
//         "Reading SDF files from a FIFO stream is not supported.");
//   }
// };
//
// SDFReader::~SDFReader() {};
//
// void SDFReader::InitializeScan(ClientContext &context,
//                                SDFReaderScanState &state) {
//   state.finished = false;
//   if (!state.file_handle || state.file_handle->path != file_handle->path) {
//     auto flags = FileFlags::FILE_FLAGS_READ;
//     state.file_handle = fs.OpenFile(file_handle->path, flags);
//   }
// }
//
// void SDFReader::Scan(SDFReaderScanState &state, DataChunk &result) {
//   while (ScanInternal(state, result)) {
//     if (result.size() > 0) {
//       break;
//     }
//     result.Reset();
//   }
// }
//
// bool SDFReader::ScanInternal(SDFReaderScanState &state, DataChunk &output) {
//   std::cout << "Scanning SDF file" << "\n";
//   if (state.finished) {
//     return false;
//   }
//   return true;
// }
//
// } // namespace duckdb
