#include <functional>
#include <iostream>
#include <vector>

template <typename element>
class disjoint_set_forest
{
public:
  struct node
  {
    node() : parent(element()), value(element()), size(0)
    {
    }
    node(const element& v) : parent(v), value(v), size(1)
    {
    }
    node(const node& n) : parent(n.parent), value(n.value), size(n.size)
    {
    }
    node& operator=(node& rhs)
    {
      return rhs;
    }
    bool operator==(const node& rhs) const
    {
      return value == rhs.value;
    }

    element parent; // Store this rather than a pointer or referece to the actual node, as later resizing will invalidate that
    element value;
    uint64_t size; // Number of descendants, including self
  };

  disjoint_set_forest(std::function<size_t(const element&)> index_map, size_t num_elements, std::function<void(const element&)> print_element)
      : _index_map(index_map), _print_element(print_element), _num_elements(num_elements)
  {
    _forest.resize(num_elements);
  }

  // IMPORTANT: v must not be in the forest.
  void make_set(const element& v)
  {
    node& n = _forest[_index_map(v)];
    n.value = v;
    n.size = 1;
    n.parent = v;
  }

  // e must be in forest
  element& find(const element& e)
  {
    node& n = _forest[_index_map(e)];

    return find(n)->value;
  }

  void merge(const element& e1, const element& e2)
  {
    node* n1 = &_forest[_index_map(e1)];
    node* n2 = &_forest[_index_map(e2)];

    n1 = find(n1);
    n2 = find(n2);

    if (n1 == n2)
    {
      return;
    }

    if (n1->size < n2->size)
    {
      n1->parent = n2->value;
      n2->size += n1->size;

      if (n2->size > _largest_tree)
      {
        _largest_tree = n2->size;
      }
    }
    else
    {
      n2->parent = n1->value;
      n1->size += n2->size;

      if (n1->size > _largest_tree)
      {
        _largest_tree = n1->size;
      }
    }
  }

  size_t get_largest_tree() const
  {
    return _largest_tree;
  }

private:
  node* find(node* n)
  {
    while (n->parent != n->value)
    {
      n->parent = get_node(n->parent)->parent;
      n = get_node(n->parent);
    }
    return n;
  }

  inline node* get_node(const element& e)
  {
    return &_forest[_index_map(e)];
  }

  std::vector<node> _forest;
  size_t _num_elements;
  std::function<size_t(const element&)> _index_map; // Must map elements to a unique index in the range [0, num_elements)
  std::function<void(const element&)> _print_element;
  size_t _largest_tree = 0;
};

/*
Problem with references approach: inplace_vector might work, but cannot be resized.
Resizing would lead to invalidating all references.
*/