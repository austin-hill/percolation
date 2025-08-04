#include <math.h>
#include <random>
#include <stdint.h>

class prng
{
public:
  prng() : _state(static_cast<uint64_t>(std::random_device{}()) << 32 | std::random_device{}())
  {
  }

  void seed()
  {
    _state = static_cast<uint64_t>(std::random_device{}()) << 32 | std::random_device{}();
  }

  inline uint64_t next_xorshift_64()
  {
    _state ^= _state << 18;
    _state ^= _state >> 31;
    _state ^= _state << 11;
    return _state;
  }

  inline uint64_t next_xorshift_64s()
  {
    _state ^= _state << 12;
    _state ^= _state >> 25;
    _state ^= _state << 27;
    return _state * 2685821657736338717;
  }

  inline uint64_t next_lcg_64()
  {
    _state = _state * 2862933555777941757 + 1;
    return (_state << 27) ^ _state;
  }

private:
  uint64_t _state;
};
