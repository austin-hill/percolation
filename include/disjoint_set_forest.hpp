#pragma once

#define force_inline inline __attribute__((always_inline))

#include <cmath>
#include <functional>
#include <iostream>
#include <print>
#include <vector>

#include "memory_mapped_vector.h"

/*
Store nodes by indices in array (64 bit), with function to retrieve value on output.
Same for parents
Store size as int, for clusters up to 2 billion in size (sign is to determine termination status - negative if cluster is on boundary).
Thus usage is 12 bytes per node.

Usage was 4*3*2+8 = 32 bytes per node.

Improved usage is now 8+4 = 12 bytes per node.
*/

template <typename element>
class disjoint_set_forest
{
public:
  /*
  Force tight packing to ensure 12 bytes per struct, otherwise will be 16. In theory tight struct packing can cause a perfomance
  penalty when accesses cross cache lines (usually 64 bytes) in a CPU. On the other hand, tight packing can improve performance by
  reducing the number of cache misses, due to more data fitting on a cache line. In testing, tight packing here seems to improve
  performance marginally, although not enough to be certain of that. Either way the memory saving is huge.
  */
  struct [[gnu::packed]] node
  {
    node() : parent_index(0), size(0)
    {
    }
    node(const element& e) : parent_index(get_index(e)), size(1)
    {
    }
    node(const node& n) : parent_index(n.parent_index), size(n.size)
    {
    }
    node& operator=(node& rhs)
    {
      return rhs;
    }

    /*
    Need this to be able to sort by size in a std::map
    NOTE: second part is to ensure different clusters are distinct, as map uses < for equality comparisons.
    */
    bool operator<(const node& rhs) const
    {
      return std::abs(size) < std::abs(rhs.size) ||
             (std::abs(size) == std::abs(rhs.size) && parent_index < rhs.parent_index); // This is eating up 20% of time
    }

    size_t parent_index; // Node index is not necessary to store - get_index is used to retrieve this
    int size;            // Number of descendants, including self
  };

  disjoint_set_forest(size_t num_elements) : _num_elements(num_elements)
  {
    _forest.resize(num_elements);
  }

  virtual size_t get_index(const element& node) const = 0; // Must map elements to a unique index in the range [0, num_elements)
  virtual element get_element(size_t index) const = 0;     // Inverse map of the above map

  virtual bool on_boundary(const element& node) const = 0;

  // IMPORTANT: element must not be in the forest.
  force_inline void make_set(const element& e)
  {
    node& n = _forest[get_index(e)];
    n.size = 2 * on_boundary(e) - 1; // Negative size implies cluster hits boundary
    n.parent_index = get_index(e);
  }

  // Element must be in forest
  force_inline element find(const element& e)
  {
    node& n = _forest[get_index(e)];

    return get_element(get_index(n));
  }

  force_inline virtual void merge(const element& e1, const element& e2)
  {
    node* n1 = &_forest[get_index(e1)];
    node* n2 = &_forest[get_index(e2)];

    n1 = find(n1);
    n2 = find(n2);

    if (n1 == n2)
    {
      return;
    }

    if (std::abs(n1->size) < std::abs(n2->size))
    {
      n1->parent_index = get_index(n2);
      // n2->size = n2->size * (1 - 2 * (n2->size > 0 && n1->size < 0)) + n1->size * (1 - 2 * (n1->size > 0 && n2->size < 0));
      // Branchless way to add sizes and set result as negative if either cluster hits the boundary (i.e. has negative size)
      n2->size = (std::abs(n1->size) + std::abs(n2->size)) * (1 - 2 * (n1->size < 0 || n2->size < 0));
    }
    else
    {
      n2->parent_index = get_index(n1);
      n1->size = (std::abs(n1->size) + std::abs(n2->size)) * (1 - 2 * (n1->size < 0 || n2->size < 0));
    }
  }

protected:
  // Find root with path halving
  force_inline node* find(node* n)
  {
    while (n->parent_index != get_index(n))
    {
      n->parent_index = _forest[n->parent_index].parent_index;
      n = &_forest[n->parent_index];
    }
    return n;
  }

  // Find root without modifying path
  force_inline const node* find_const(const node* n) const
  {
    while (n->parent_index != get_index(n))
    {
      n = &_forest[n->parent_index];
    }
    return n;
  }

  force_inline node* get_node(const element& e)
  {
    return &_forest[get_index(e)];
  }

  force_inline size_t get_index(const node* n) const
  {
    return n - &_forest[0];
  }

  std::vector<node> _forest;
  size_t _num_elements;
};

/*
Problem with references approach: inplace_vector might work, but cannot be resized.
Resizing would lead to invalidating all references.
*/