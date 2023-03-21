/**
 * @file initializer.h
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief 
 * @version 0.0.1
 * @date 2023-03-11
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef QUADRATUBE_MODEL_INITIALIZER_H_
#define QUADRATUBE_MODEL_INITIALIZER_H_

#include <string>

#include "core/math.h"
#include "model/system.h"
#include "metadata.h"

namespace ModelInitializer {

typedef struct {
  int m;       ///< parastichy number m
  int n;       ///< parastichy number n
  int repeat;  ///< repeat times
  int direction; ///< direction of dislocations
  int glide;   ///< steps of glide
  int climb;   ///< steps of climb
  int bn;      ///< begin position
  double rest_len;  ///< rest length
} Parameters;

class Initializer {
  public:
    inline Initializer(ModelSystem& system): __system(system) {}
    void init(Parameters init_parameter);

  private:
    /// @brief flags of function
    enum ObjectType {
      Bond1, Bond2, Curvatrue
    };

    /// @brief to surrport reuse of functions
    inline CoreMath::View<CoreMath::Array<int>>& node_adjacents(ObjectType tp) {
      if (tp == Bond1)
        return __system.node_adjacents_bonds1_;
      if (tp == Bond2)
        return __system.node_adjacents_bonds2_;
      return __system.node_adjacents_curvature_;
    }

    /// @brief make sure don't input Curvatrue
    inline CoreMath::View<CoreMath::Pair<int>>& bond_relations(ObjectType tp) {
      if (tp == Bond1)
        return __system.bond_relations1_;
      return __system.bond_relations2_;
    }

    /// @brief calculate index from 2-dimensional coordinates
    inline int flat(int i, int j) {
      j += __para.m * (i/__para.n);
      i %= __para.n;
      if (i < 0) {
        i += __para.n;
        j -= __para.m;
      }
      j %= __para.m * __para.repeat;
      if (j < 0)
        j += __para.m * __para.repeat;
      return i + __para.n*j;
    }

    /// @brief find nodes connected by bond(a, b) oppsite to ck
    int dual(ObjectType tp, int ck, int a, int b);
    inline int rot(ObjectType tp, int ck, int a, int b) {
      auto i =  node_adjacents(tp)[a].find(b);
      if (i ==  node_adjacents(tp)[a].end())
        return -1;
      auto il = (i == node_adjacents(tp)[a].begin()) ?
          node_adjacents(tp)[a].end()-1 : i-1;
      auto ir = (i == node_adjacents(tp)[a].end()-1) ?
          node_adjacents(tp)[a].begin() : i+1;
      return (*il == ck) ? *ir : *il;
    }

    /// @brief remove bond
    inline void remove(ObjectType tp, CoreMath::Pair<int> bond) {
      node_adjacents(tp)[bond[0]].erase(
          node_adjacents(tp)[bond[0]].find(bond[1]));
      node_adjacents(tp)[bond[1]].erase(
          node_adjacents(tp)[bond[1]].find(bond[0]));
      if (tp != Curvatrue)
        bond_relations(tp).remove(bond);
    }

    /// @brief remove node node and related bonds
    void remove(int node);
    /// @brief bond, n1 is insert position of bond[0], n2 is insert position of bond[1]
    inline void insert(ObjectType tp, CoreMath::Pair<int> bond, int* n1, int* n2) {
      if (tp != Curvatrue)
        bond_relations(tp).push_back(bond);
      node_adjacents(tp)[bond[0]].insert(n1, bond[1]);
      node_adjacents(tp)[bond[1]].insert(n2, bond[0]);
    }
    
    /// @brief add p to last of node_positions_, return it's position
    inline int insert(CoreMath::Vector p) {
      __system.node_positions_.push_back(p);
      return __system.node_positions_.size() - 1;
    }
    
    /// @brief between bond(i, j) and bond(i, k), position to insert
    inline int* between(ObjectType tp, int i, int j, int k) {
      int* m = node_adjacents(tp)[i].find(j);
      int* n = node_adjacents(tp)[i].find(k);
      return (m == node_adjacents(tp)[i].begin()) ?
          n+1 : ((n == node_adjacents(tp)[i].begin()) ?
              m+1 : ((n>m) ? n : m));
    }

    Parameters __para;
    ModelSystem& __system;
};

} // namespace ModelInitializer

#endif // QUADRATUBE_MODEL_INITIALIZER_H_
