/**
 * Implementation of treap.hpp
 */

#include "treap.hpp"
#include <stdio.h>

Treap::Treap::Treap() {
  this->root = nullptr;
}

Treap::Treap::Treap(Node *root) {
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
  return nullptr;
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
    } else {
      return n;
    }
  }
  return nullptr;
}





