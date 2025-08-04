#include <chrono>
#include <iostream>
#include <stdint.h>

// add ifndef

class timer
{
public:
  timer() : _time_ns(0)
  {
  }

  inline void start()
  {
    _start = std::chrono::high_resolution_clock::now();
  }
  inline void restart()
  {
    _time_ns = 0;
    _start = std::chrono::high_resolution_clock::now();
  }
  inline void stop()
  {
    _time_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - _start).count();
  }
  uint64_t get_ns()
  {
    return _time_ns;
  }
  uint64_t get_ms()
  {
    return _time_ns / 1000;
  }
  void print_ns()
  {
    std::cout << "Time elapsed: " << _time_ns << "ns" << std::endl;
  }
  void print_us()
  {
    std::cout << "Time elapsed: " << _time_ns / 1e3 << "ms" << std::endl;
  }
  void print_ms()
  {
    std::cout << "Time elapsed: " << _time_ns / 1e6 << "ms" << std::endl;
  }

private:
  std::chrono::high_resolution_clock::time_point _start;
  uint64_t _time_ns;
};
