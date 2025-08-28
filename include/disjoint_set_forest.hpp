#pragma once

#include <functional>
#include <iostream>
#include <print>
#include <vector>

/*
Store nodes by indices in array (64 bit), with function to retrieve value on output.
Same for parents
Store size as uint32, for clusters up to 4 billion in size.
Thus usage is 12 bytes per node.

Usage was 4*3*2+8 = 32 bytes per node.

Improved usage is now 8+4 = 12 bytes per node.
*/

template <typename element>
class disjoint_set_forest
{
public:
  struct node
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
      return size < rhs.size || (size == rhs.size && parent_index < rhs.parent_index);
    }

    size_t parent_index; // Node index is not necessary to store
    uint32_t size;       // Number of descendants, including self
  };

  disjoint_set_forest(size_t num_elements) : _num_elements(num_elements)
  {
    _forest.resize(num_elements);
  }

  virtual size_t get_index(const element& node) = 0; // Must map elements to a unique index in the range [0, num_elements)
  virtual element get_element(size_t index) = 0;     // Inverse map of the above map

  // IMPORTANT: v must not be in the forest.
  void make_set(const element& e)
  {
    node& n = _forest[get_index(e)];
    n.size = 1;
    n.parent_index = get_index(e);
  }

  // e must be in forest
  element& find(const element& e)
  {
    node& n = _forest[get_index(e)];

    return get_element(get_index(n));
  }

  virtual void merge(const element& e1, const element& e2)
  {
    node* n1 = &_forest[get_index(e1)];
    node* n2 = &_forest[get_index(e2)];

    n1 = find(n1);
    n2 = find(n2);

    if (n1 == n2)
    {
      return;
    }

    if (n1->size < n2->size)
    {
      n1->parent_index = get_index(n2);
      n2->size += n1->size;
    }
    else
    {
      n2->parent_index = get_index(n1);
      n1->size += n2->size;
    }
  }

protected:
  node* find(node* n)
  {
    while (n->parent_index != get_index(n))
    {
      n->parent_index = _forest[n->parent_index].parent_index;
      n = &_forest[n->parent_index];
    }
    return n;
  }

  const node* find_const(const node* n)
  {
    while (n->parent_index != get_index(n))
    {
      n = &_forest[n->parent_index];
    }
    return n;
  }

  node* get_node(const element& e)
  {
    return &_forest[get_index(e)];
  }

  size_t get_index(const node* n) const
  {
    return n - &*_forest.begin();
  }

  std::vector<node> _forest;
  size_t _num_elements;
};

/*
Problem with references approach: inplace_vector might work, but cannot be resized.
Resizing would lead to invalidating all references.
*/