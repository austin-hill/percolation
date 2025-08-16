#pragma once

#include <algorithm>
#include <boost/container_hash/hash.hpp>
#include <limits>
#include <print>
#include <queue>
#include <random>
#include <set>
#include <stdint.h>

#include "gnuplot-iostream.h"

#include "pcg_extras.hpp"
#include "pcg_random.hpp"
#include "percolation.h"
#include "timer.h"

template <size_t cube_size>
class cubic_site_percolation
{
  /*
  Can possibly use a disjoint set union data structure.
  Loop over all nodes in cube.
  Draw edges to previous nodes (three different possible edges): If connecting any two, merge them.
  Expected complexity: Hopefully in most cases we don't have to do much merging. Either way, should be amortized O(n*ackerman^-1(n)).
  */
public:
  static size_t get_index(const std::tuple<int, int, int>& node)
  {
    return std::get<0>(node) + std::get<1>(node) * cube_size + std::get<2>(node) * cube_size * cube_size;
  }

  static void print_node(const std::tuple<int, int, int>& node)
  {
    std::cout << "Node: (" << std::get<0>(node) << ", " << std::get<1>(node) << ", " << std::get<2>(node) << ")" << std::endl;
  }

  cubic_site_percolation(double p)
      : _p(p), _bound(std::numeric_limits<uint64_t>::max() * p), _nodes(get_index, cube_size * cube_size * cube_size, print_node),
        _rng(pcg_extras::seed_seq_from<std::random_device>{})
  {
    _gp << "set xrange [0:" << cube_size << "]" << std::endl;
    _gp << "set yrange [0:" << cube_size << "]" << std::endl;
    _gp << "set zrange [0:" << cube_size << "]" << std::endl;

    // _gp << "splot NaN" << std::endl;
  }

  inline std::vector<std::tuple<int, int, int>> get_previous(const std::tuple<int, int, int>& node)
  {
    return {
        {    std::get<0>(node),     std::get<1>(node), std::get<2>(node) - 1},
        {    std::get<0>(node), std::get<1>(node) - 1,     std::get<2>(node)},
        {std::get<0>(node) - 1,     std::get<1>(node),     std::get<2>(node)},
    };
  }

  bool generate_cluster()
  {
    std::cout << "Generating cluster..." << std::endl;

    for (int i = 0; i < cube_size; ++i)
    {
      for (int j = 0; j < cube_size; ++j)
      {
        for (int k = 0; k < cube_size; ++k)
        {
          // Make set and merge with previous (up to three) if there is an edge between them
          std::tuple<int, int, int> new_node(i, j, k);
          _nodes.make_set(new_node);

          for (const auto& node : get_previous(new_node))
          {
            if (std::get<0>(node) < 0 || std::get<1>(node) < 0 || std::get<2>(node) < 0)
            {
              continue;
            }
            if (_rng() < _bound)
            {
              _nodes.merge(node, new_node);
            }
          }
        }
      }
    }

    std::cout << "Largest cluster: " << _nodes.get_largest_tree() << std::endl;

    auto largest_graph = _nodes.get_largest_graph_nodes();

    _gp << "set title 'test'\n";
    _gp << "splot NaN" << std::endl;

    const auto largest_clusters = _nodes.get_clusters_sorted(10000);
    const std::vector<std::string> colours = {"red", "orange", "yellow", "green", "blue", "purple", "pink", "brown", "grey"};
    size_t colour_index = 0;
    for (auto it = largest_clusters.crbegin(); it != largest_clusters.crend() && colour_index < 10; ++it, ++colour_index)
    { //  rgb " << colours[colour_index] << "
      _gp << "replot" << _gp.file1d(it->second) << "u 1:2:3 with points pt 7 ps 0.001 title 'test'" << std::endl;
      sleep(1);
    }
    std::println("Number of clusters of size at least 10000: {}", largest_clusters.size());

    // _gp << "splot" << _gp.file1d(largest_graph) << "u 1:2:3 with points pt 7 ps 0.001 title 'test'" << std::endl;

    return true;
  }

private:
  const double _p;
  const uint64_t _bound;
  percolation<std::tuple<int, int, int>> _nodes; // All of the nodes in the cluster visited so far
  pcg64_fast _rng;
  Gnuplot _gp;
  // timer _tm;
};