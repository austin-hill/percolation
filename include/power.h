#pragma once

#include <stdint.h>
#include <type_traits>

template <typename int_base, std::enable_if_t<std::is_integral<int_base>::value, bool> = true, typename int_exp,
          std::enable_if_t<std::is_unsigned<int_exp>::value, bool> = true>
int_base ipow(int_base base, int_exp exp)
{
  int_base result = 1;
  for (;;)
  {
    if (exp & 1)
      result *= base;
    exp >>= 1;
    if (!exp)
      break;
    base *= base;
  }

  return result;
}

template <uint64_t base, uint8_t exp>
struct ipow_tmp
{
  static constexpr uint64_t value = base * ipow_tmp<base, exp - 1>::value;
};

template <uint64_t base>
struct ipow_tmp<base, 0>
{
  static constexpr uint64_t value = 1;
};
