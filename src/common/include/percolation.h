#pragma once

#include <algorithm>
#include <limits>
#include <map>
#include <print>
#include <queue>
#include <stdint.h>

#include "gnuplot-iostream.h"

#include "disjoint_set_forest.hpp"
#include "timer.h"

// In here, have definitions of merge etc which keep track of largest clusters. Have functions which can get clusters etc.
// Add enable if too.
// Do we want largest N clusters, or every cluster of size > N? Probably the latter?

template <typename element>
class percolation : public disjoint_set_forest<element>
{
  using node = typename percolation<element>::node;

  /*
  Can possibly use a disjoint set union data structure.
  Loop over all nodes in cube.
  Draw edges to previous nodes (three different possible edges): If connecting any two, merge them.
  Expected complexity: Hopefully in most cases we don't have to do much merging. Either way, should be amortized O(n*ackerman^-1(n)).
  */
public:
  struct root_node
  {
    root_node(const element& e, uint32_t size, bool terminated) : value(e), size(size), terminated(terminated)
    {
    }

    element value;
    uint32_t size;
    bool terminated;
  };

  percolation(size_t num_elements) : disjoint_set_forest<element>(num_elements)
  {
  }

  virtual bool on_boundary(const element& node) = 0;

  virtual void print_node(const std::tuple<int, int, int>& node) = 0;

  /* std::vector<std::vector<element>> get_clusters(size_t minimum_size)
  {
    std::vector<std::vector<element>> clusters;
    std::unordered_map<size_t, size_t> root_index_to_cluster_index_map;

    for (auto& n : this->_forest)
    {
      const auto& root_node = *this->find(&n);
      if (root_node.size > minimum_size)
      {
        const auto root_index = this->_index_map(root_node.value);
        if (root_index_to_cluster_index_map.contains(root_index))
        {
          clusters[root_index_to_cluster_index_map.at(root_index)].push_back(n.value);
        }
        else
        {
          root_index_to_cluster_index_map.emplace(root_index, clusters.size());
          clusters.push_back(std::vector<element>({n.value}));
        }
      }
    }

    return clusters;
  } */

  std::map<node, std::pair<std::vector<element>, bool>> get_clusters_sorted(size_t minimum_size)
  {
    std::map<node, std::pair<std::vector<element>, bool>> clusters;

    for (size_t index = 0; index < this->_forest.size(); ++index)
    {
      const node& root = *this->find(&this->_forest[index]);

      if (root.size > minimum_size)
      {
        if (clusters.contains(root))
        {
          clusters.at(root).first.push_back(this->get_element(index));
          clusters.at(root).second |= on_boundary(this->get_element(index));
        }
        else
        {
          clusters.emplace(
              root, std::pair<std::vector<element>, bool>{std::vector<element>({this->get_element(index)}), on_boundary(this->get_element(index))});
        }
      }
    }

    return clusters;
  }
};