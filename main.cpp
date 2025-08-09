/* Aim:
 calculate percolation phase transition value as accurately as possible.
 Idea:
 Use PCG prng (or faster one if needed).
 Use multiple threads to grow clusters in parallel.
 Start from below and increase.
 How do we grow clusters?
 How do we figure out termination?
 All we need to know is all of the paths OUT of a given cube which are connected
 to the cluster. Problem: how do we revisit cubes?

 To start, just do single thread version, no cubes, etc.

 What about union find? Make loads of small clusters in separate cubes and merge
 them?
*/

#include <algorithm>
#include <boost/container_hash/hash.hpp>
#include <limits>
#include <queue>
#include <random>
#include <set>
#include <stdint.h>
#include <unordered_set>

#include "gnuplot-iostream.h"

#include "include/disjoint_set_forest.hpp"
#include "include/flat_hash_map.hpp"
#include "include/pcg_extras.hpp"
#include "include/pcg_random.hpp"
#include "include/timer.h"

/* void print_node(const std::tuple<int, int, int>& node)
{
  std::cout << "Node: (" << std::get<0>(node) << ", " << std::get<1>(node) << ", " << std::get<2>(node) << ")" << std::endl;
} */

class percolation
{
public:
  percolation(double p) : _p(p), _bound(std::numeric_limits<uint64_t>::max() * p), _rng(pcg_extras::seed_seq_from<std::random_device>{})
  {
    _neighbours.push({0, 0, 0});

    _gp << "set xrange [-50:50]" << std::endl;
    _gp << "set yrange [-50:50]" << std::endl;
    _gp << "set zrange [-50:50]" << std::endl;

    _gp << "splot NaN" << std::endl;
  }

  void print_node(const std::tuple<int, int, int>& node)
  {
    std::cout << "Node: (" << std::get<0>(node) << ", " << std::get<1>(node) << ", " << std::get<2>(node) << ")"
              << " Neighbours size: " << _neighbours.size() << std::endl;
  }

  inline std::vector<std::tuple<int, int, int>> get_possible_neighbours(const std::tuple<int, int, int>& node)
  {
    return {
        {    std::get<0>(node),     std::get<1>(node), std::get<2>(node) + 1},
        {    std::get<0>(node),     std::get<1>(node), std::get<2>(node) - 1},
        {    std::get<0>(node), std::get<1>(node) + 1,     std::get<2>(node)},
        {    std::get<0>(node), std::get<1>(node) - 1,     std::get<2>(node)},
        {std::get<0>(node) + 1,     std::get<1>(node),     std::get<2>(node)},
        {std::get<0>(node) - 1,     std::get<1>(node),     std::get<2>(node)},
    };
  }

  bool generate_cluster()
  {
    std::cout << "Generating cluster..." << std::endl;
    // Almost all time is spent checking and inserting into _nodes.
    while (_neighbours.size() != 0)
    {
      // Add next node to queue
      auto next_node = _neighbours.front();
      _neighbours.pop();
      _nodes.insert(next_node);

      // print_node(next_node);

      // Generate its neighbours
      for (auto n : get_possible_neighbours(next_node))
      {
        bool res = _nodes.count(n);
        if (!res && _rng() < _bound)
        {
          _neighbours.push(n);
          _gp << "set arrow from " << std::get<0>(next_node) << ", " << std::get<1>(next_node) << ", " << std::get<2>(next_node) << " to "
              << std::get<0>(n) << ", " << std::get<1>(n) << ", " << std::get<2>(n) << "nohead lw 0.2" << std::endl;
        }
      }

      if (_nodes.size() % 1000 == 0)
      {
        _gp << "refresh" << std::endl;
      }

      if (_nodes.size() > 1e5)
      {
        std::cout << "Failed to terminate: " << _neighbours.size() << std::endl;

        /* std::vector<std::tuple<int, int, int>> nb;
        while (_neighbours.size() > 0)
        {
          nb.push_back(_neighbours.front());
          _neighbours.pop();
        } */

        /* Gnuplot gp;

        gp << "set title 'test'\n";

        gp << "splot" << gp.file1d(nb) << "u 1:2:3 with points pt 7 ps 0.1 title 'test'" << std::endl; */

        //_tm.print_ms();
        return false;
      }
    }
    std::cout << "Generated cluster of size " << _nodes.size() << std::endl;
    return true;
  }

private:
  const double _p;
  const uint64_t _bound;
  ska::flat_hash_set<std::tuple<int, int, int>, boost::hash<std::tuple<int, int, int>>> _nodes; // All of the nodes in the cluster visited so far
  std::queue<std::tuple<int, int, int>> _neighbours; // All of the immediate neighbours of the so far visited cluster
  pcg64_fast _rng;
  Gnuplot _gp;
  // timer _tm;
};

template <size_t cube_size>
class percolation1
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

  percolation1(double p)
      : _p(p), _bound(std::numeric_limits<uint64_t>::max() * p), _nodes(get_index, cube_size * cube_size * cube_size, print_node),
        _rng(pcg_extras::seed_seq_from<std::random_device>{})
  {
    _gp << "set xrange [0:" << cube_size << "]" << std::endl;
    _gp << "set yrange [0:" << cube_size << "]" << std::endl;
    _gp << "set zrange [0:" << cube_size << "]" << std::endl;

    _gp << "splot NaN" << std::endl;
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
    _gp << "splot" << _gp.file1d(largest_graph) << "u 1:2:3 with points pt 7 ps 0.01 title 'test'" << std::endl;

    return true;
  }

private:
  const double _p;
  const uint64_t _bound;
  disjoint_set_forest<std::tuple<int, int, int>> _nodes; // All of the nodes in the cluster visited so far
  pcg64_fast _rng;
  Gnuplot _gp;
  // timer _tm;
};

int main()
{
  percolation1<400> p(0.248); // p(0.2488125);
  // percolation p(0.248);

  timer tm;
  tm.start();
  p.generate_cluster();
  tm.stop();
  tm.print_ms();

  return 0;
}

/*
Memory usage:
4 bytes per int *3 plus 8 bytes size plus 1 byte pointer, so total of 21 bytes per node.

So 500*500*500 probably uses about 2.7GB of RAM

100*100*100 uses about 21MB of RAM.
*/