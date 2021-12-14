/**
 * @file treap.hpp
 * @brief Treap implementation.
 *
 * A Treap (Tree + Heap) is a binary tree where each node
 * contains one element x with key(x) and prio(x). Treaps
 * have binary serach tree and heap properties.
 *
 * @author Atamert Rahma, atamertiel@gmail.com
 */

#ifndef TREAP_HPP
#define TREAP_HPP

#include <limits>
#include <cmath>

/**
 * Tools for Treaps
 */
namespace Treap {

  struct Pair {
    /**
     * Should be given unique. (We assume this!)
     * Should be non-negative if the node is not an empty leaf.
     */
    int key;
    //! Non-negative unique Priority.
    int prio;
  };

  /**
   * @struct Node
   * Represents a treap node. An empty leaf is represented
   * when isEmpty is set.
   */
  struct Node {
    Pair p = {-1, std::numeric_limits<int>::min()};
    Node *left = nullptr;
    Node *right = nullptr;
    Node *parent = nullptr;
    bool isEmpty = true;
  };

  /**
   * Some inline functions for better readability in the implementation.
   */
  inline int key(Node *n) {
    return n->p.key;
  }
  inline int prio(Node *n) {
    return n->p.prio;
  }
  inline Node* parent(Node *n) {
    return n->parent;
  }
  inline Node* left(Node *n) {
    return n->left;
  }
  inline Node* right(Node *n) {
    return n->right;
  }
  inline bool isLeaf(Node *n) {
    return !(right(n) || left(n)) && !n->isEmpty;
  }
  inline bool isEmptyLeaf(Node *n) {
    return n->isEmpty
        && key(n) == -1
        && prio(n) == std::numeric_limits<int>::min();
  }
  /**
   * Adds left and right empty leaves to the given node.
   */
  inline void addEmptyLeaves(Node *n) {
    /* left empty leaf */
    Node *l = new Node();
    n->left = l;
    l->parent = n;
    /* right empty leaf */
    Node *r = new Node();
    n->right = r;
    r->parent = n;
  }

  /**
   * @class Treap
   *
   * A Treap (Tree + Heap) is a binary tree where each node
   * contains one element x with key(x) and prio(x) where
   * the following properties hold:
   *
   * 1) Search Tree property: For each element x:
   *   - elements y in the left subtree of x satisfy: key(y) < key(x)
   *   - elements y in the right subtree of x satisfy: key(y) > key(x)
   *   - we assume that the key values are pairwise distinct.
   *
   * 2) Heap property: For all elements x,y:
   *   - if y is a child of x, then prio(y) > prio(x)
   *   - all priorities must be pairwise distinct.
   *
   * For elements x_1, .. , x_n with key(x_i) and prio(x_i), there exists
   * a unique Treap. And it has the structure that would result if
   * elements were inserted in the order of their priorities.
   */
  class Treap {
  public:
    Treap();
    Treap(Node *root = nullptr);

    /* Iterative implementation */
    Node* searchKey(int k);
    /* Recursive implementation */
    Node* searchKey(Node *n, int k);

    void insertNode(Node *n);
    void deleteNode(Node *n);
    Node* min();
    Node* max();
    Node* suc(unsigned int k);
    Node* pre(unsigned int k);

  private:
    Node *root = nullptr;
    void rotateLeft(Node *n);
    void rotateRight(Node *n);
  };

  Node* min(Treap *t);
  Node* max(Treap *t);
  Node* list(Treap *t);
  void retUnion(Treap *t1, Treap *t2);
  void retSplit(Treap *t, unsigned int k, Treap *t1, Treap *t2);
}

#endif
