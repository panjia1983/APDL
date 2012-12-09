#ifndef CSPACE_LEARNING_LSHFOREST_H
#define CSPACE_LEARNING_LSHFOREST_H

#include <map>

namespace cspace_learning
{

template<typename LSH, typename Key>
class ForestIndex2
{
public:
  typedef typename LSH::Parameter Parameter;
  typedef typename LSH::Domain Domain;

protected:

  struct Tree
  {
    std::vector<LSH> lshs; // the hash functions

    struct Node
    {
      size_t size; // total # points in subtree
      std::vector<Node*> children;
      std::vector<Key> data;

      Node() : size(0) {}

      ~Node()
      {
        for(unsigned int i = 0; i < children.size(); ++i)
        {
          delete children[i];
        }
      }

      bool empty() const
      {
        return size == 0;
      }

      template<typename Accessor>
      void insert(Tree* tree, unsigned int depth, Key key, Accessor& accessor)
      {
        ++size;
        if(children.empty())
        {
          data.push_back(key);
          if(depth < tree->lshs.size() && data.size() > 1)
          {
            LSH& lsh = tree->lshs[depth];
            if(lsh.getRange() == 0)
              throw std::logic_error("LSH with unlimited hash value cannot be used in the forest LSH.");
            children.resize(lsh.getRange());
            for(unsigned int i = 0; i < data.size(); ++i)
            {
              Key key = data[i];
              unsigned int h = lsh(accessor(key));
              if(children[h] == 0)
                children[h] = new Node();
              children[h]->insert(tree, depth+1, key, accessor);
            }
            data.clear();
          }
        }
        else
        {
          unsigned int h = tree->lshs[depth](accessor(key));
          if(children[h] == 0)
            children[h] = new Node();
          children[h]->insert(tree, depth+1, key, accessor);
        }
      }

      template<typename Scanner>
      void scan(Scanner& scanner) const
      {
        if(!children.empty())
        {
          for(unsigned int i = 0; i < children.size(); ++i)
          {
            Node* n = children[i];
            if(n) n->scan(scanner);
          }
        }

        if(!data.empty())
        {
          for(unsigned int i = 0; i < data.size(); ++i)
          {
            scanner(data[i]);
          }
        }
      }
    };

    Node* root;

  public:
    Tree() : root(0) {}

    template<typename RNG>
    void reset(const Parameter& param, RNG& rng, unsigned int depth)
    {
      lshs.resize(depth);
      for(unsigned int i = 0; i < lshs.size(); ++i)
      {
        lshs[i].reset(param, rng);
      }

      root = new Node();
    }

    ~Tree()
    {
      if(root) delete root;
    }

    template<typename Accessor>
    void insert(Key key, Accessor& accessor)
    {
      root->insert(this, 0, key, accessor);
    }

    void lookup(Domain val, std::vector<const Node*>& nodes) const
    {
      const Node* cur = root;
      unsigned int depth = 0;
      nodes.clear();
      for(;;)
      {
        nodes.push_back(cur);
        if(cur->children.empty()) break;
        unsigned int h = lshs[depth](val);
        cur = cur->children[h];
        if(!cur) break;
        ++depth;
      }
    }
  };


  std::vector<Tree> trees;

public:
  ForestIndex2() {}

  template<typename RNG>
  void init(const Parameter& param, RNG& rng, unsigned int L, unsigned int depth)
  {
    trees.resize(L);
    for(unsigned int i = 0; i < trees.size(); ++i)
      trees[i].reset(param, rng, depth);
  }

  template<typename Accessor>
  void insert(Key key, Accessor& accessor)
  {
    for(unsigned int i = 0; i < trees.size(); ++i)
    {
      trees[i].insert(key, accessor);
    }
  }

  template<typename Scanner>
  void query(unsigned int M, Scanner& scanner) const
  {
    std::vector<std::vector<const typename Tree::Node *> > list(trees.size());
    for(unsigned int i = 0; i < trees.size(); ++i)
    {
      trees[i].lookup(scanner.query(), list[i]);
    }

    // find the minimal depth convering at least M points
    unsigned int d = 0;
    for(;;)
    {
      unsigned int s = 0;
      for(unsigned int i = 0; i < list.size(); ++i)
      {
        if(d < list[i].size())
          s += list[i][d]->size;
      }

      if(s < M) break;
      ++d;
    }

    if(d > 0) --d;

    // recursively scan the nodes
    for(unsigned int i = 0; i < list.size(); ++i)
    {
      if(d < list[i].size())
        list[i][d]->scan(scanner);
    }
  }
};


template<typename LSH, typename Key>
class ForestIndex
{
public:
  typedef typename LSH::Parameter Parameter;
  typedef typename LSH::Domain Domain;

protected:

  struct Tree
  {
    std::vector<LSH> lshs; // the hash functions

    struct Node
    {
      typedef std::map<unsigned int, Node*> Children;

      size_t size;
      Children children;
      std::vector<Key> data;

      Node() : size(0) {}

      ~Node()
      {
        for(typename Children::iterator i = children.begin(); i != children.end(); ++i)
          delete i->second;
      }

      bool empty() const
      {
        return size == 0;
      }

      void insert(Tree* tree, unsigned int depth, Key key, Domain value)
      {
        ++size;
        if(children.empty())
        {
          data.push_back(key);
          if(depth < tree->lshs.size() && data.size() > 1)
          {
            LSH& lsh = tree->lshs[depth];
            for(unsigned int i = 0; i < data.size(); ++i)
            {
              Key key = data[i];
              unsigned int h = lsh(value);
              if(children.find(h) == children.end())
                children[h] = new Node();
              children[h]->insert(tree, depth+1, key, value);
            }
            data.clear();
          }
        }
        else
        {
          unsigned int h = tree->lshs[depth](value);
          if(children.find(h) == children.end())
            children[h] = new Node();
          children[h]->insert(tree, depth+1, key, value);
        }
      }

      template<typename Scanner>
      void scan(Scanner& scanner) const
      {
        if(!children.empty())
        {
          for(typename Children::const_iterator i = children.begin(); i != children.end(); ++i)
          {
            Node* n = i->second;
            n->scan(scanner);
          }
        }

        if(!data.empty())
        {
          for(unsigned int i = 0; i < data.size(); ++i)
          {
            scanner(data[i]);
          }
        }
      }
    };

    Node* root;

  public:
    Tree() : root(0) {}

    template<typename RNG>
    void reset(const Parameter& param, RNG& rng, unsigned int depth)
    {
      lshs.resize(depth);
      for(unsigned int i = 0; i < lshs.size(); ++i)
      {
        lshs[i].reset(param, rng);
      }

      root = new Node();
    }

    ~Tree()
    {
      if(root) delete root;
    }

    void insert(Key key, Domain value)
    {
      root->insert(this, 0, key, value);
    }

    void lookup(Domain val, std::vector<const Node*>& nodes) const
    {
      Node* cur = root;
      unsigned int depth = 0;
      nodes.clear();
      for(;;)
      {
        nodes.push_back(cur);
        if(cur->children.empty()) break;
        unsigned int h = lshs[depth](val);
        if(cur->children.find(h) == cur->children.end())
        {
          cur = 0;
          break;
        }
        else
          cur = cur->children[h];
        ++depth;
      }
    }
  };


  std::vector<Tree> trees;

public:
  ForestIndex() {}

  template<typename RNG>
  void init(const Parameter& param, RNG& rng, unsigned int L, unsigned int depth)
  {
    trees.resize(L);
    for(unsigned int i = 0; i < trees.size(); ++i)
      trees[i].reset(param, rng, depth);
  }

  void insert(Key key, Domain value)
  {
    for(unsigned int i = 0; i < trees.size(); ++i)
    {
      trees[i].insert(key, value);
    }
  }

  template<typename Scanner>
  void query(unsigned int M, Scanner& scanner) const
  {
    std::vector<std::vector<const typename Tree::Node *> > list(trees.size());
    for(unsigned int i = 0; i < trees.size(); ++i)
    {
      trees[i].lookup(scanner.query(), list[i]);
    }

    // find the minimal depth convering at least M points
    unsigned int d = 0;
    for(;;)
    {
      unsigned int s = 0;
      for(unsigned int i = 0; i < list.size(); ++i)
      {
        if(d < list[i].size())
          s += list[i][d]->size;
      }

      if(s < M) break;
      ++d;
    }

    if(d > 0) --d;

    // recursively scan the nodes
    for(unsigned int i = 0; i < list.size(); ++i)
    {
      if(d < list[i].size())
        list[i][d]->scan(scanner);
    }
  }
};


}

#endif
