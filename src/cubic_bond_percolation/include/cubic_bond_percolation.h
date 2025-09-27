#pragma once

#define force_inline inline __attribute__((always_inline))

#include <set>
#include <stdint.h>

#include "gnuplot-iostream.h"

#include "pcg_extras.hpp"
#include "pcg_random.hpp"
#include "percolation.h"
#include "power.h"

// GNU plot has its limitations here. Do not waste too much time fiddling with it, will probably write something proper later anyway.

class cubic_bond_percolation : public percolation<std::tuple<int, int, int>>
{
  /*
  Can possibly use a disjoint set union data structure.
  Loop over all nodes in cube.
  Draw edges to previous nodes (three different possible edges): If connecting any two, merge them.
  Expected complexity: Hopefully in most cases we don't have to do much merging. Either way, should be amortized O(n*ackerman^-1(n)).
  */
public:
  cubic_bond_percolation(double p);

  size_t get_index(const std::tuple<int, int, int>& node) const override;
  std::tuple<int, int, int> get_element(size_t index) const override;
  bool on_boundary(const std::tuple<int, int, int>& node) const override;

  void generate_clusters();
  void generate_clusters_parallel(uint8_t max_num_threads);
  void plot_clusters(uint32_t min_cluster_size, size_t max_num_clusters = 10) const;
  void plot_central_clusters(uint32_t min_cluster_size, size_t central_cube_size = 64, size_t max_num_clusters = 10) const;

private:
  void generate_merge_clusters_recursive(uint8_t max_num_threads, int start_i, int end_i);
  void generate_clusters_parallel_thread(int start_i, int end_i);
  void merge_clusters_slices(int i);

  const double _p;
  static constexpr uint8_t _cube_pow = 10;
  static constexpr uint32_t _cube_size = ipow_tmp<2, _cube_pow>::value;
  const uint64_t _bound;
  pcg64_fast _rng;
  mutable Gnuplot _gp;
};

// Need to speed this up... Maybe write in assembly by hand
force_inline size_t cubic_bond_percolation::get_index(const std::tuple<int, int, int>& node) const
{
  return static_cast<size_t>(std::get<0>(node)) | (static_cast<size_t>(std::get<1>(node) << _cube_pow)) |
         (static_cast<size_t>(std::get<2>(node)) << (2 * _cube_pow));
}

force_inline std::tuple<int, int, int> cubic_bond_percolation::get_element(size_t index) const
{
  std::tuple<int, int, int> element;
  std::get<2>(element) = index >> (2 * _cube_pow);
  std::get<1>(element) = (index >> _cube_pow) - (std::get<2>(element) << _cube_pow);
  std::get<0>(element) = index - (std::get<1>(element) << _cube_pow) - (std::get<2>(element) << (2 * _cube_pow));
  return element;
}

force_inline bool cubic_bond_percolation::on_boundary(const std::tuple<int, int, int>& node) const
{
  return std::get<0>(node) == 0 || std::get<0>(node) == _cube_size - 1 || std::get<1>(node) == 0 || std::get<1>(node) == _cube_size - 1 ||
         std::get<2>(node) == 0 || std::get<2>(node) == _cube_size - 1;
}
