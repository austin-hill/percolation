#pragma once

#include <set>
#include <stdint.h>

#include "gnuplot-iostream.h"

#include "pcg_extras.hpp"
#include "pcg_random.hpp"
#include "percolation.h"

// GNU plot has its limitations here. Do not waste too much time fiddling with it, will probably write something proper later anyway.

class cubic_site_percolation : percolation<std::tuple<int, int, int>>
{
  /*
  Can possibly use a disjoint set union data structure.
  Loop over all nodes in cube.
  Draw edges to previous nodes (three different possible edges): If connecting any two, merge them.
  Expected complexity: Hopefully in most cases we don't have to do much merging. Either way, should be amortized O(n*ackerman^-1(n)).
  */
public:
  cubic_site_percolation(double p, uint8_t cube_size_pow);

  size_t get_index(const std::tuple<int, int, int>& node) override;
  std::tuple<int, int, int> get_element(size_t index) override;
  bool on_boundary(const std::tuple<int, int, int>& node) override;

  void print_node(const std::tuple<int, int, int>& node) override
  {
    std::cout << "Node: (" << std::get<0>(node) << ", " << std::get<1>(node) << ", " << std::get<2>(node) << ")" << std::endl;
  }

  inline std::vector<std::tuple<int, int, int>> get_previous(const std::tuple<int, int, int>& node);

  bool generate_clusters();
  bool generate_clusters_parallel(uint8_t max_num_threads);
  void plot_clusters(uint32_t min_cluster_size, size_t max_num_clusters = 10);

private:
  bool generate_merge_clusters_recursive(uint8_t max_num_threads, int start_i, int end_i);
  bool generate_clusters_parallel_thread(int start_i, int end_i);
  bool merge_clusters_slices(int i);

  const double _p;
  const uint32_t _cube_size;
  const uint8_t _cube_pow;
  const uint64_t _bound;
  pcg64_fast _rng;
  Gnuplot _gp;
};