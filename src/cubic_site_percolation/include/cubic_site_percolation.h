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

#include "colour_names.h"
#include "pcg_extras.hpp"
#include "pcg_random.hpp"
#include "percolation.h"
#include "timer.h"

// GNU plot has its limitations here. Do not waste too much time fiddling with it, will probably write something proper later anyway.

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

  static bool on_boundary(const std::tuple<int, int, int>& node)
  {
    return std::get<0>(node) == 0 || std::get<0>(node) == cube_size - 1 || std::get<1>(node) == 0 || std::get<1>(node) == cube_size - 1 ||
           std::get<2>(node) == 0 || std::get<2>(node) == cube_size - 1;
  }

  static void print_node(const std::tuple<int, int, int>& node)
  {
    std::cout << "Node: (" << std::get<0>(node) << ", " << std::get<1>(node) << ", " << std::get<2>(node) << ")" << std::endl;
  }

  cubic_site_percolation(double p)
      : _p(p), _bound(std::numeric_limits<uint64_t>::max() * p), _nodes(get_index, cube_size * cube_size * cube_size, print_node, on_boundary),
        _rng(pcg_extras::seed_seq_from<std::random_device>{})
  {
    _gp << "set xrange [0:" << cube_size << "]" << std::endl;
    _gp << "set yrange [0:" << cube_size << "]" << std::endl;
    _gp << "set zrange [0:" << cube_size << "]" << std::endl;
    _gp << "set view equal xyz" << std::endl; // Force grid to be square
    _gp << "unset border" << std::endl;
    _gp << "unset xtics" << std::endl;
    _gp << "unset ytics" << std::endl;
    _gp << "unset ztics" << std::endl;
    _gp << "set key outside right top samplen 2 spacing .7 font ',8' tc rgb 'grey40'" << std::endl;
  }

  inline std::vector<std::tuple<int, int, int>> get_previous(const std::tuple<int, int, int>& node)
  {
    return {
        {    std::get<0>(node),     std::get<1>(node), std::get<2>(node) - 1},
        {    std::get<0>(node), std::get<1>(node) - 1,     std::get<2>(node)},
        {std::get<0>(node) - 1,     std::get<1>(node),     std::get<2>(node)},
    };
  }

  bool generate_clusters()
  {
    std::cout << "Generating clusters..." << std::endl;

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

    return true;
  }

  void plot_clusters(uint32_t min_cluster_size, size_t max_num_clusters = 10)
  {
    max_num_clusters = std::min(colour_names.size(), max_num_clusters);
    _gp << "set title tc rgb 'grey40' 'Percolation, p=" << std::setprecision(8) << _p << " Cube size=" << cube_size << "'" << std::endl;

    const auto largest_clusters = _nodes.get_clusters_sorted(min_cluster_size);
    size_t count = 0;

    std::println("Number of clusters of size at least {}: {}", min_cluster_size, largest_clusters.size());

    for (auto it = largest_clusters.crbegin(); it != largest_clusters.crend() && count < max_num_clusters; ++it, ++count)
    {
      _gp << ((count == 0) ? "splot" : "replot") << _gp.file1d(it->second.first) << "u 1:2:3:(0.03) with points lc rgb '" << colour_names[count]
          << "' pt 7 ps variable title 'Cluster " << count + 1 << " (" << it->second.first.size() << " points)"
          << ((it->second.second) ? "(still growing)" : "(terminated)") << "'" << ((count == max_num_clusters - 1) ? "; pause mouse close" : "")
          << std::endl;
    }
  }

private:
  const double _p;
  const uint64_t _bound;
  percolation<std::tuple<int, int, int>> _nodes; // All of the nodes in the cluster visited so far
  pcg64_fast _rng;
  Gnuplot _gp;
  // timer _tm;
};