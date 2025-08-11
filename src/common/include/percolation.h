#pragma once

#include <algorithm>
#include <limits>
#include <queue>
#include <set>
#include <stdint.h>

#include "gnuplot-iostream.h"

#include "disjoint_set_forest.hpp"
#include "timer.h"

// In here, have definitions of merge etc which keep track of largest clusters. Have functions which can get clusters etc.
// Add enable if too.

template <typename element>
class percolation : public disjoint_set_forest<element>
{
  using node = typename disjoint_set_forest<element>::node;

  /*
  Can possibly use a disjoint set union data structure.
  Loop over all nodes in cube.
  Draw edges to previous nodes (three different possible edges): If connecting any two, merge them.
  Expected complexity: Hopefully in most cases we don't have to do much merging. Either way, should be amortized O(n*ackerman^-1(n)).
  */
public:
  percolation(std::function<size_t(const element&)> index_map, size_t num_elements, std::function<void(const element&)> print_element)
      : disjoint_set_forest<element>(index_map, num_elements, print_element)
  {
  }

  void merge(const element& e1, const element& e2) override
  {
    node* n1 = &this->_forest[this->_index_map(e1)];
    node* n2 = &this->_forest[this->_index_map(e2)];

    n1 = this->find(n1);
    n2 = this->find(n2);

    if (n1 == n2)
    {
      return;
    }

    if (n1->size < n2->size)
    {
      n1->parent = n2->value;
      n2->size += n1->size;

      if (n2->size > this->_forest[_largest_tree_root].size)
      {
        _largest_tree_root = this->_index_map(n2->value);
      }
    }
    else
    {
      n2->parent = n1->value;
      n1->size += n2->size;

      if (n1->size > this->_forest[_largest_tree_root].size)
      {
        _largest_tree_root = this->_index_map(n1->value);
      }
    }
  }

  size_t get_largest_tree() const
  {
    return this->_forest[_largest_tree_root].size;
  }

  std::vector<element> get_largest_graph_nodes()
  {
    std::vector<element> largest_graph;
    const node& root_node = this->_forest[_largest_tree_root];

    for (auto& n : this->_forest)
    {
      if (*this->find(&n) == root_node)
      {
        largest_graph.push_back(n.value);
      }
    }

    return largest_graph;
  }

private:
  size_t _largest_tree_root = 0;
};