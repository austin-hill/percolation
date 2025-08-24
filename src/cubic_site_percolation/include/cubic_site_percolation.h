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
#include "power.h"
#include "timer.h"

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
  size_t get_index(const std::tuple<int, int, int>& node) override
  {
    return std::get<0>(node) + (std::get<1>(node) << _cube_pow) + (std::get<2>(node) << (2 * _cube_pow));
  }

  std::tuple<int, int, int> get_element(size_t index) override
  {
    std::tuple<int, int, int> element;
    std::get<2>(element) = index >> (2 * _cube_pow);
    std::get<1>(element) = (index >> _cube_pow) - (std::get<2>(element) << _cube_pow);
    std::get<0>(element) = index - (std::get<1>(element) << _cube_pow) - (std::get<2>(element) << (2 * _cube_pow));
    return element;
  }

  bool on_boundary(const std::tuple<int, int, int>& node) override
  {
    return std::get<0>(node) == 0 || std::get<0>(node) == _cube_size - 1 || std::get<1>(node) == 0 || std::get<1>(node) == _cube_size - 1 ||
           std::get<2>(node) == 0 || std::get<2>(node) == _cube_size - 1;
  }

  void print_node(const std::tuple<int, int, int>& node) override
  {
    std::cout << "Node: (" << std::get<0>(node) << ", " << std::get<1>(node) << ", " << std::get<2>(node) << ")" << std::endl;
  }

  cubic_site_percolation(double p, uint8_t cube_size_pow)
      : _p(p), _cube_size(ipow(static_cast<uint32_t>(2), cube_size_pow)), _cube_pow(cube_size_pow), _bound(std::numeric_limits<uint64_t>::max() * p),
        percolation(ipow(static_cast<uint32_t>(2), static_cast<uint8_t>(cube_size_pow * 3))), _rng(pcg_extras::seed_seq_from<std::random_device>{})
  {
    _gp << "set xrange [0:" << _cube_size << "]" << std::endl;
    _gp << "set yrange [0:" << _cube_size << "]" << std::endl;
    _gp << "set zrange [0:" << _cube_size << "]" << std::endl;
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

    for (int i = 0; i < _cube_size; ++i)
    {
      for (int j = 0; j < _cube_size; ++j)
      {
        for (int k = 0; k < _cube_size; ++k)
        {
          // Make set and merge with previous (up to three) if there is an edge between them
          std::tuple<int, int, int> new_node(i, j, k);
          this->make_set(new_node);

          for (const auto& node : get_previous(new_node))
          {
            if (std::get<0>(node) < 0 || std::get<1>(node) < 0 || std::get<2>(node) < 0)
            {
              continue;
            }
            if (_rng() < _bound)
            {
              this->merge(node, new_node);
            }
          }
        }
      }
    }

    std::println("Done");

    return true;
  }

  void plot_clusters(uint32_t min_cluster_size, size_t max_num_clusters = 10)
  {
    max_num_clusters = std::min(colour_names.size(), max_num_clusters);
    _gp << "set title tc rgb 'grey40' 'Percolation, p=" << std::setprecision(8) << _p << " Cube size=" << _cube_size << "'" << std::endl;

    const auto largest_clusters = this->get_clusters_sorted(min_cluster_size);
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
  const uint32_t _cube_size;
  const uint8_t _cube_pow;
  const uint64_t _bound;
  // percolation<std::tuple<int, int, int>> _nodes; // All of the nodes in the cluster visited so far
  pcg64_fast _rng;
  Gnuplot _gp;
  // timer _tm;
};