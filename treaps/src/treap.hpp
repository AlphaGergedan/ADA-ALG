#ifndef TREAP_HPP
#define TREAP_HPP

namespace Treap {

  struct Pair {
    /* Unique Key */
    unsigned int key;
    /* Unique Priority */
    unsigned int prio;
  };

  struct Node {
    Pair p;
    Node *left;
    Node *right;
  };

  inline unsigned int key(Node *n) {
    return n->p.key;
  }

  inline unsigned int prio(Node *n) {
    return n->p.prio;
  }

  class Treap {
  public:
    Treap();
    Treap(Node *root = nullptr);
    void print();

    Node* searchKey(unsigned int k);
    Node* searchKey(Node *n, unsigned int k);

    void insertNode(Node *n);
    void deleteNode(Node *n);
    Node* min();
    Node* max();
    Node* suc(unsigned int k);
    Node* pre(unsigned int k);

    int getMaxDepth();
    int getMinDepth();
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
