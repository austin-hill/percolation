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