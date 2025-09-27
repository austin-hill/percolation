#include <cstring>
#include <format>
#include <stdexcept>
#include <string>

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

template <typename T>
class memory_mapped_vector
{
public:
  memory_mapped_vector() : _size(_page_size) // TODO: this could be unnecessarily large depending on T
  {
    _fd = open("/tmp/", O_RDWR | O_TMPFILE, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (_fd == -1)
    {
      throw std::runtime_error("Failed to create temporary file");
    }

    if (ftruncate(_fd, _size * sizeof(T)) == -1)
    {
      throw std::runtime_error("Failed resize temporary file");
    }

    _data = static_cast<T*>(mmap(nullptr, _size * sizeof(T), PROT_READ | PROT_WRITE, MAP_PRIVATE, _fd, 0));
    if (_data == MAP_FAILED)
    {
      throw std::runtime_error(std::format("Failed to map file: Code {}", errno));
    }
  }

  ~memory_mapped_vector()
  {
    if (munmap((void*)_data, _size * sizeof(T)) == -1)
    {
      close(_fd);
      std::println("Failed to unmap temporary file");
    }

    close(_fd);
  }

  void resize(size_t new_size)
  {
    if (new_size > _size)
    {
      if (ftruncate(_fd, new_size * sizeof(T)) == -1)
      {
        throw std::runtime_error("Failed resize temporary file");
      }

      _data = static_cast<T*>(mremap(_data, _size * sizeof(T), new_size * sizeof(T), MREMAP_MAYMOVE));

      if (_data == MAP_FAILED)
      {
        throw std::runtime_error(std::format("Failed to remap file: Code {}", errno));
      }

      _size = new_size;
    }
  }

  T& operator[](size_t n) noexcept
  {
    assert(n < _size);
    return _data[n];
  }

  const T& operator[](size_t n) const noexcept
  {
    assert(n < _size);
    return _data[n];
  }

  constexpr size_t size() const noexcept
  {
    return _size;
  }

private:
  T* _data;
  size_t _size;

  int _fd;

  static constexpr size_t _page_size = 4096;
};
