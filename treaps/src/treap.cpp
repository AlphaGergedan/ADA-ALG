/**
 * @file treap.cpp
 * @brief Implementation of treap.hpp
 * @author Atamert Rahma, atamertiel@gmail.com
 */

#include "treap.hpp"
#include <iostream>

/**
 * Initializes an empty leaf.
 */
Treap::Treap::Treap() {
  this->root = new Node();
}

/**
 * Initializes a Treap from the given node (root).
 * @param root    Pointer to a node, must not be empty or have parent.
 */
Treap::Treap::Treap(Node *root) {
  if (!root) {
    /* Initializes an empty leaf. */
    this->root = new Node();
  }
  else {
    if (parent(root)) {
      std::cerr << "The given node has a parent node so it cannot be a root.";
    }
    else if (isEmptyLeaf(root)) {
      std::cerr << "The given node is an empty leaf!";
    }
    else {
      this->root = root;
    }
  }
}

/**
 * Search for the node that has the given key value.
 *
 * Time: O(#elements-on-the-search-path)
 * Expected Running Time: O(log n)
 *
 * @param k     Integer key value to search in the treap.
 * @return v    Pointer to the Node with the key value k or
 *              the empty leaf where the search ends
 *              (fails because the element does not exist).
 */
Treap::Node* Treap::Treap::searchKey(int k) {
  Node *v = this->root;
  while (!isEmptyLeaf(v)) {
    /* element found */
    if (key(v) == k) {
      return v;
    }
    /* element resides in the right subtree of v */
    else if (key(v) < k) {
      v = v->right;
    }
    /* element resides in the left subtree of v */
    else {
      v = v->left;
    }
  }
  /* element not found. Return empty leaf. */
  return v;
}

/**
 * Search for the node that has the given key value.
 * Recursive implementation.
 *
 * Time: O(#elements-on-the-search-path)
 * Expected Running Time: O(log n)
 *
 * @param n     Pointer to the search root
 *              (where the search starts).
 * @param k     Integer key value to search in the treap.
 * @return n    Pointer to the Node with the key value k or
 *              the empty leaf where the search ends
 *              (fails because the element does not exist).
 */
Treap::Node* Treap::Treap::searchKey(Node *n, int k) {
  if (!isEmptyLeaf(n)) {
    /* element resides in the right subtree of n */
    if (k > key(n)) {
      return searchKey(n->right, k);
    }
    /* element resides in the left subtree of n */
    else if (k < key(n)) {
      return searchKey(n->left, k);
    }
    /* element found */
    else {
      return n;
    }
  }
  /* element not found. Return empty leaf. */
  return n;
}

void Treap::Treap::insertNode(Node *n) {
  /* first search for the element */
  Node *x = this->searchKey(key(n));

  /* x must be a leaf node */
  if (!x) {
    /* Insert the node as a leaf. This may only violate the heap property. */
    x = n;

    /* fix the heap property */
    while (prio(x->parent) ) {
      // TODO
    }
  } else {
    std::cout << "The node is already present in the treap!" << std::endl;
  }
}





