// #include "sdf_scanner/buffered_sdf_reader.hpp"
// #include "duckdb/common/assert.hpp"
// #include "duckdb/common/file_opener.hpp"
// #include "duckdb/common/unique_ptr.hpp"
// #include <utility>
// namespace duckdb {
//
// SDFBufferHandle::SDFBufferHandle(idx_t buffer_index_p, idx_t readers_p,
//                                  AllocatedData &&buffer_p, idx_t
//                                  buffer_size_p)
//     : buffer_index(buffer_index_p), readers(readers_p),
//       buffer(std::move(buffer_p)), buffer_size(buffer_size_p) {}
//
// BufferedSDFReader::BufferedSDFReader(ClientContext &context,
//                                      BufferedSDFReaderOptions options_p,
//                                      string file_name_p)
//     : context(context), options(std::move(options_p)),
//       file_name(std::move(file_name_p)), buffer_index(0), thrown(false) {}
//
// void BufferedSDFReader::OpenSDFFile() {
//   lock_guard<mutex> guard(lock);
//   if (!IsOpen()) {
//     auto &fs = FileSystem::GetFileSystem(context);
//     auto regular_file_handle =
//         fs.OpenFile(file_name, FileFlags::FILE_FLAGS_READ);
//     file_handle = make_uniq<SDFFileHandle>(std::move(regular_file_handle),
//                                            BufferAllocator::Get(context));
//   }
//   Reset();
// }
//
// bool SDFFileHandle::Read(char *pointer, idx_t &read_size, idx_t
// requested_size,
//                          bool &file_done, bool sample_run) {
//
//   // TODO: Implement SDF reading code here instead
//   // D_ASSERT(requested_size != 0);
//   // if (last_read_requested) {
//   //   return false;
//   // }
//   //
//   // if (can_seek) {
//   //   read_size = ReadInternal(pointer, requested_size);
//   //   read_position += read_size;
//   // } else if (sample_run) { // Cache the buffer
//   //   read_size = ReadInternal(pointer, requested_size);
//   //   if (read_size > 0) {
//   //     cached_buffers.emplace_back(allocator.Allocate(read_size));
//   //     memcpy(cached_buffers.back().get(), pointer, read_size);
//   //   }
//   //   cached_size += read_size;
//   //   read_position += read_size;
//   // } else {
//   //   read_size = 0;
//   //   if (!cached_buffers.empty() || read_position < cached_size) {
//   //     read_size += ReadFromCache(pointer, requested_size, read_position);
//   //   }
//   //   if (requested_size != 0) {
//   //     read_size += ReadInternal(pointer, requested_size);
//   //   }
//   // }
//   //
//   // if (read_size == 0) {
//   //   last_read_requested = true;
//   //   file_done = true;
//   // }
//   std::cout << "Reading SDF File..." << std::endl;
//   return true;
// }
//
//
// } // namespace duckdb
