#pragma once
#include "duckdb/main/client_context.hpp"
namespace duckdb {

struct SDFFileHandle {
public:
  SDFFileHandle(unique_ptr<FileHandle> file_handle, Allocator &allocator);

  bool IsOpen() const;
  void Close();

  void Reset();
  bool RequestedReadsComplete();
  bool LastReadRequested() const;

  idx_t FileSize() const;
  idx_t Remaining() const;

  bool CanSeek() const;
  bool IsPipe() const;

  FileHandle &GetHandle();

  //! The next two functions return whether the read was successful
  bool GetPositionAndSize(idx_t &position, idx_t &size, idx_t requested_size);
  bool Read(char *pointer, idx_t &read_size, idx_t requested_size,
            bool &file_done, bool sample_run);
  //! Read at position optionally allows passing a custom handle to read from,
  //! otherwise the default one is used
  void ReadAtPosition(char *pointer, idx_t size, idx_t position,
                      bool &file_done, bool sample_run,
                      optional_ptr<FileHandle> override_handle = nullptr);

private:
  idx_t ReadInternal(char *pointer, const idx_t requested_size);
  // idx_t ReadFromCache(char *&pointer, idx_t &size, idx_t &position);

private:
  //! The SDF file handle
  unique_ptr<FileHandle> file_handle;
  Allocator &allocator;

  //! File properties
  const bool can_seek;
  const idx_t file_size;

  //! Read properties
  idx_t read_position;
  atomic<idx_t> requested_reads;
  atomic<idx_t> actual_reads;
  atomic<bool> last_read_requested;

  //! Cached buffers for resetting when reading stream
  vector<AllocatedData> cached_buffers;
  idx_t cached_size;
};

class PlainSDFReader {
public:
  PlainSDFReader(ClientContext &context, string file_name);
  void OpenSDFFile();
  void Reset();

  bool HasFileHandle() const;
  bool IsOpen() const;
  const string &GetFileName() const;
  SDFFileHandle &GetFileHandle() const;

  mutable mutex lock;

private:
  ClientContext &context;

  //! File name
  const string file_name;
  //! File handle
  unique_ptr<SDFFileHandle> file_handle;
};
} // namespace duckdb
