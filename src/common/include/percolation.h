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

template <typename element>
class percolation : public disjoint_set_forest<element>
{
  /*
  Can possibly use a disjoint set union data structure.
  Loop over all nodes in cube.
  Draw edges to previous nodes (three different possible edges): If connecting any two, merge them.
  Expected complexity: Should be amortized O(n*ackerman^-1(n)).
  */
public:
  using node = typename percolation<element>::node;

  percolation(size_t num_elements) : disjoint_set_forest<element>(num_elements)
  {
  }

  std::map<node, std::vector<element>> get_clusters_sorted(size_t minimum_size) const
  {
    std::map<node, std::vector<element>> clusters;

    for (size_t index = 0; index < this->_forest.size(); ++index)
    {
      const node& root = *this->find_const(&this->_forest[index]);

      if (std::abs(root.size) >= minimum_size)
      {
        if (clusters.contains(root))
        {
          clusters.at(root).push_back(this->get_element(index));
        }
        else
        {
          clusters.emplace(root, std::vector<element>({this->get_element(index)}));
        }
      }
    }

    return clusters;
  }
};