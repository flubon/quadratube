/**
 * @file math.h
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief Array and other mathematical algrithoms
 * @version 0.0.1
 * @date 2023-03-07
 * 
 * @copyright Copyright (c) 2023
 */
#ifndef QUADRATUBE_CORE_MATH_H_
#define QUADRATUBE_CORE_MATH_H_

#include <Kokkos_Vector.hpp>
#include <Kokkos_Random.hpp>

#define PI 3.141592653589793
//! Boltzmann constant
#define K_B 1.380649e-23

namespace CoreMath {

/**
 * @class Array
 * @brief Array type of capacity N
 * @details All members satisfy definition of vector in C++ standard library,
 *     except 'find()'. Be careful because we never check bounds of array due to 
 *     performance.
 * 
 * @tparam T 
 */
template <typename T, size_t N = 8>
class Array {
  public:
    /// @brief Constructor
    inline Array(size_t n = 0): __len(n), __data{T()} {}

    /// @brief Iterator
    typedef T* iterator;
    inline iterator begin() { return &__data[0]; }
    inline iterator end() { return &__data[__len]; }

    /// @brief Capacity
    inline size_t size() { return __len; }
    inline void resize(size_t n) { __len = n; }

    /// @brief Element access
    inline T& operator[](int i) { return __data[i]; }
    inline T* data() { return __data; }
    
    /// @brief Modifiers
    inline void push_back(T in) { __data[__len++] = in;}
    inline T pop_back() { return __data[--__len];}

    iterator insert(iterator pos, const T& val) {
      for (iterator i = __data[__len]; i != pos; i--)
        *i = *(i - 1);
      *pos = val;
      __len++;
      return pos;
    }

    iterator erase(iterator pos) {
      for (iterator i = pos + 1; i != __data[__len]; i++)
        *(i - 1) = *i;
      __len--;
      return pos;
    }

    /// @brief Compare operators
    bool operator==(Array<T> p) {
      if (__len != p.__len)
        return false;
      for (int i=0; i<__len; i++)
        if (__data[i] != __data[i])
          return false;
      return true;
    }

    /// @brief Additional functions
    iterator find(T from) { 
      for (iterator i = __data; i != &__data[__len]; i++)
        if (*i == from)
          return i;
      return &__data[__len];
    }

  private:
    size_t __len;
    T __data[N];
}; // class Array


/**
 * @class View
 * @brief Dynamic array
 * @details We called view instead of vector here because vector is 3 dimensional
 *     vector here.
 * 
 * @tparam T 
 */
template <typename T>
class View : Kokkos::vector<T> {
  
}; // class View


/**
 * @class Vector
 * @brief 3 dimensional vector
 * @details The following functions are operations of vectors.
 */
class Vector {
  public:
    inline Vector(): __data{0} {}
    inline Vector(double i, double j, double k): __data{i, j, k} {}
    inline double& operator[](int i) { return __data[i]; }

    /// @brief add, subtract and self-add
    inline Vector operator+(Vector p) {
      return Vector(__data[0] + p[0], __data[1] + p[1], __data[2] + p[2]);
    }
    inline Vector operator-() { return Vector(-__data[0], -__data[1], -__data[2]); }
    inline Vector operator-(Vector p) {
      return Vector(__data[0] - p[0], __data[1] - p[1], __data[2] - p[2]);
    }
    inline void operator+=(Vector p) {
      __data[0] += p[0]; __data[1] += p[1]; __data[2] += p[2];
    }

    /// @brief quantity product
    inline Vector operator*(double n) {
      return Vector(n * __data[0], n * __data[1], n * __data[2]);
    }
    
    /// @brief 1/n times vector
    inline Vector operator/(double n) {
      return Vector(__data[0] / n, __data[1] / n, __data[2] / n);
    }

    /// @brief dot product of two vectors
    inline double operator*(Vector p) {
      return __data[0] * p[0] + __data[1] * p[1] + __data[2] * p[2];
    }

    /// @brief compare two vectors
    inline bool operator==(Vector p) {
      return __data[0]==p[0] && __data[1]==p[1] && __data[2]==p[2];
    }

  private:
    double __data[3];
}; // class Vector

/// @brief quantity product for n on the left hand
inline Vector operator*(double n, Vector p) { return p * n; }

/// @brief cross product of two vectors
inline Vector cross(Vector p1, Vector p2) {
  return Vector( p1[1]*p2[2] - p1[2]*p2[1], p1[2]*p2[0] - p1[0]*p2[2],
    p1[0]*p2[1] - p1[1]*p2[0]);
}

/// @brief length of vector p
inline double mod(Vector p) { return Kokkos::sqrt(p * p); }


/**
 * @class Pair
 * @brief A pair of something
 * @details A pair of nodes is a bond. We can have a pair of dislocations too.
 * 
 * @tparam T 
 */
template <typename T>
class Pair {
  public:
    inline Pair(): __data{0} {}
    inline Pair(T a, T b): __data{a, b} {}
    inline T& operator[](int i) { return __data[i]; }
    
    /// @brief If node a is connected by this bond, return the other node. Return -1 else.
    inline T find(T a) {
      return (__data[0] == a) ? __data[1] : ((__data[1] == a) ? __data[0] : -1);
    }

    /// @brief replace a with b，return false if failed，otherwise true.
    inline bool replace(T a, T b) {
      if (__data[0] == a) {
        __data[0] = b;
        return true;
      }
      if (__data[1] == a) {
        __data[1] = b;
        return true;
      }
      return false;
    }

    /// @brief compare
    inline bool operator==(Pair<T> b) {
      return (__data[0] == b[0] && __data[1] == b[1]) || (__data[1] == b[0] && __data[0] == b[1]);
    }

    /// @brief add and substrct, for 2 dimensional coordinates
    inline Pair<T> operator-(Pair<T> b) {
      return Pair<T>(__data[0] - b[0], __data[1] - b[1]);
    }
    inline Pair<T> operator+(Pair<T> b) {
      return Pair<T>(__data[0] + b[0], __data[1] + b[1]);
    }

  private:
    T __data[2];
}; // class Pair


/**
 * @brief Random pool for vector generation
 * 
 */
class Pool {
  public:
    inline Pool() : __pool(Kokkos::Timer().seconds()) {}
    inline Vector gen_vector(double end) {
      auto generator = __pool.get_state();
      double r = generator.drand(0, end);
      double phi = generator.drand(0, 2*PI);
      double theta = generator.drand(-PI/2, PI/2);
      return Vector(r * Kokkos::cos(phi) * Kokkos::sin(theta),
        r * Kokkos::sin(phi) * Kokkos::sin(theta), r * Kokkos::cos(theta));
    }

  private:
    Kokkos::Random_XorShift64_Pool<> __pool;
}; // class Pool

} // namespace CoreMath

#endif // QUADRATUBE_CORE_MATH_H_
