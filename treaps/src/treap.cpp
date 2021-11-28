/**
 * Implementation of treap.hpp
 */

#include "treap.hpp"
#include <iostream>

Treap::Treap::Treap() {
  this->root = new Node();
  this->root->left = nullptr;
  this->root->right = nullptr;
  this->root->parent= nullptr;
}

Treap::Treap::Treap(Node *root) {
  if (parent(root)) {
    std::cerr << "The given node has a parent node so it cannot be a root.";
  }
  this->root = root;
}

void Treap::Treap::print() {
  int depth = this->getMaxDepth();
  // TODO
}

/**
 *
 * Time: O(#elements-on-the-search-path)
 */
Treap::Node* Treap::Treap::searchKey(unsigned int k) {
  Node *v = this->root;
  while (v) {
    if (key(v) == k) {
      return v;
    }
    else if (key(v) < k) {
      /* search tree property */
      v = v->right;
    }
    else {
      v = v->left;
    }
  }
  /* element not found */
  return v;
}

/**
 * Recursive search implementation.
 *
 * Time: O(#elements-on-the-search-path)
 *
 * @param n     Pointer to the search root
 * @param k     key that we are searching
 * @return      Pointer to the node that contains the key 'k' or nullptr
 */
Treap::Node* Treap::Treap::searchKey(Node *n, unsigned int k) {
  if (n) {
    if (k > key(n)) {
      searchKey(n->right, k);
    }
    else if (k < key(n)) {
      searchKey(n->left, k);
    }
    else {
      return n;
    }
  }
  return nullptr;
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





