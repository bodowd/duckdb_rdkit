#include "sdf_scanner/plain_sdf_reader.hpp"
#include "duckdb/common/file_system.hpp"
#include <iostream>
namespace duckdb {

SDFFileHandle::SDFFileHandle(unique_ptr<FileHandle> file_handle_p,
                             Allocator &allocator_p)
    : file_handle(std::move(file_handle_p)), allocator(allocator_p),
      can_seek(file_handle->CanSeek()), file_size(file_handle->GetFileSize()),
      read_position(0), requested_reads(0), actual_reads(0),
      last_read_requested(false), cached_size(0) {}

bool SDFFileHandle::IsOpen() const { return file_handle != nullptr; }

void SDFFileHandle::Close() {
  if (IsOpen() && !file_handle->IsPipe()) {
    file_handle->Close();
    file_handle = nullptr;
  }
}

void SDFFileHandle::Reset() {
  D_ASSERT(RequestedReadsComplete());
  read_position = 0;
  requested_reads = 0;
  actual_reads = 0;
  last_read_requested = false;
  if (IsOpen() && CanSeek()) {
    file_handle->Reset();
  }
}

bool SDFFileHandle::RequestedReadsComplete() {
  return requested_reads == actual_reads;
}

bool SDFFileHandle::LastReadRequested() const { return last_read_requested; }

idx_t SDFFileHandle::FileSize() const { return file_size; }

idx_t SDFFileHandle::Remaining() const { return file_size - read_position; }

bool SDFFileHandle::CanSeek() const { return can_seek; }

bool SDFFileHandle::IsPipe() const { return file_handle->IsPipe(); }

FileHandle &SDFFileHandle::GetHandle() { return *file_handle; }

bool SDFFileHandle::GetPositionAndSize(idx_t &position, idx_t &size,
                                       idx_t requested_size) {
  D_ASSERT(requested_size != 0);
  if (last_read_requested) {
    return false;
  }

  position = read_position;
  size = MinValue<idx_t>(requested_size, Remaining());
  read_position += size;

  requested_reads++;
  if (size == 0) {
    last_read_requested = true;
  }

  return true;
}

bool SDFFileHandle::Read(char *pointer, idx_t &read_size, idx_t requested_size,
                         bool &file_done, bool sample_run) {
  std::cout << "Read from SDFFileHandle" << std::endl;
  return true;
}

void SDFFileHandle::ReadAtPosition(char *pointer, idx_t size, idx_t position,
                                   bool &file_done, bool sample_run,
                                   optional_ptr<FileHandle> override_handle) {
  std::cout << "Read at position not implemented" << std::endl;
}

idx_t SDFFileHandle::ReadInternal(char *pointer, const idx_t requested_size) {
  idx_t total_read_size = 0;
  while (total_read_size < requested_size) {
    auto read_size = file_handle->Read(pointer + total_read_size,
                                       requested_size - total_read_size);
    if (read_size == 0) {
      break;
    }
    total_read_size += read_size;
  }
  return total_read_size;
}

PlainSDFReader::PlainSDFReader(ClientContext &context, string file_name)
    : context(context), file_name(std::move(file_name)) {}

void PlainSDFReader::OpenSDFFile() {
  lock_guard<mutex> guard(lock);
  if (!IsOpen()) {
    auto &fs = FileSystem::GetFileSystem(context);
    auto regular_file_handle =
        fs.OpenFile(file_name, FileFlags::FILE_FLAGS_READ);
    file_handle = make_uniq<SDFFileHandle>(std::move(regular_file_handle),
                                           BufferAllocator::Get(context));
  }
  Reset();
}

void PlainSDFReader::Reset() {
  if (HasFileHandle()) {
    file_handle->Reset();
  }
}

bool PlainSDFReader::HasFileHandle() const { return file_handle != nullptr; }

bool PlainSDFReader::IsOpen() const {
  if (HasFileHandle()) {
    return file_handle->IsOpen();
  }
  return false;
}

const string &PlainSDFReader::GetFileName() const { return file_name; }

SDFFileHandle &PlainSDFReader::GetFileHandle() const {
  D_ASSERT(HasFileHandle());
  return *file_handle;
}

} // namespace duckdb
