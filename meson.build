project(
  'percolation',
  'cpp',
  default_options : [
    'cpp_std=c++23',
  ]
)

add_global_arguments('-O2', language : 'cpp')

library('flat_hash_map', 'external/flat_hash_map.hpp')
library('pcg', ['external/pcg_random.hpp', 'external/pcg_extras.hpp'])

library('timer', 'include/timer.h')
library('xorshift', 'include/xorshift.h')


include_dirs = include_directories([
  'external',
  'include',
  'src/common/include',
  'src/cubic_site_percolation/include',
])

executable(
  'percolation',
  'main.cpp',
  include_directories : include_dirs,
  link_args: [
    '-L/usr/lib',
    '-lboost_filesystem',
    '-lboost_system',
    '-lboost_iostreams',
  ]
)