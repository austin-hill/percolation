project(
  'percolation',
  'cpp',
  default_options: [
    'cpp_std=c++23',
  ],
)

add_global_arguments('-O2', language: 'cpp')
add_global_arguments('-flto', language: 'cpp')
add_global_arguments('-std=c++23', language: 'cpp')
add_global_arguments('-g', language: 'cpp')
add_global_arguments('-Werror', language: 'cpp')
add_global_arguments('-Wall', language: 'cpp')
add_global_arguments('-Wno-sign-compare', language: 'cpp')

flat_hash_map_lib = static_library('flat_hash_map', 'external/flat_hash_map.hpp')
flat_hash_map = declare_dependency(link_with: flat_hash_map_lib, include_directories: 'external')

pcg_lib = static_library('pcg', ['external/pcg_random.hpp', 'external/pcg_extras.hpp'])
pcg = declare_dependency(link_with: pcg_lib, include_directories: 'external')

timer_lib = static_library('timer', 'include/timer.h')
timer = declare_dependency(link_with: timer_lib, include_directories: 'include')

xorshift_lib = static_library('xorshift', 'include/xorshift.h')
xorshift = declare_dependency(link_with: xorshift_lib, include_directories: 'include')

gnuplot_link_args = [
  '-L/usr/lib',
  '-lboost_filesystem',
  '-lboost_system',
  '-lboost_iostreams',
]

executable(
  'naive_approach',
  'src/naive_approach/main.cpp',
  dependencies: [
    pcg,
    timer,
    flat_hash_map,
  ],
  link_args: gnuplot_link_args,
)

executable(
  'cubic_site_percolation',
  'src/cubic_site_percolation/cubic_site_percolation.cpp',
  include_directories: [
    'src/common/include',
    'src/cubic_site_percolation/include',
  ],
  dependencies: [
    pcg,
    timer,
  ],
  link_args: gnuplot_link_args,
)