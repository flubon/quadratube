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

#include <string>

#include <Kokkos_DualView.hpp>
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
    KOKKOS_INLINE_FUNCTION
    Array(size_t n = 0): __len(n), __data{T()} {}

    /// @brief Iterator
    using iterator = T*;
    KOKKOS_INLINE_FUNCTION
    iterator begin() { return &__data[0]; }
    KOKKOS_INLINE_FUNCTION
    iterator end() { return &__data[__len]; }

    /// @brief Capacity
    KOKKOS_INLINE_FUNCTION
    size_t size() const { return __len; }
    KOKKOS_INLINE_FUNCTION 
    void resize(size_t n) { __len = n; }

    /// @brief Element access
    KOKKOS_INLINE_FUNCTION
    T& operator[](int i) { return __data[i]; }
    KOKKOS_INLINE_FUNCTION 
    T operator[](int i) const { return __data[i]; }
    KOKKOS_INLINE_FUNCTION
    T* data() { return __data; }
    
    /// @brief Modifiers
    KOKKOS_INLINE_FUNCTION
    void push_back(T in) { __data[__len++] = in;}
    KOKKOS_INLINE_FUNCTION
    T pop_back() { return __data[--__len];}

    KOKKOS_FUNCTION
    iterator insert(iterator pos, const T& val) {
      for (iterator i = __data + __len; i > pos; i--)
        *i = *(i - 1);
      *pos = val;
      __len++;
      return pos;
    }

    KOKKOS_FUNCTION
    iterator erase(iterator pos) {
      for (iterator i = pos + 1; i < __data + __len; i++)
        *(i - 1) = *i;
      __len--;
      return pos;
    }
    
    /// @brief Additional functions
    KOKKOS_FUNCTION
    iterator find(T from) { 
      for (iterator i = __data; i < __data + __len; i++)
        if (*i == from)
          return i;
      return __data + __len;
    }

    /// @brief Compare operators
    KOKKOS_FUNCTION
    bool operator==(const Array<T>& p) const {
      if (__len != p.__len)
        return false;
      for (int i=0; i<__len; i++)
        if (__data[i] != __data[i])
          return false;
      return true;
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
class View : public Kokkos::DualView<T*> {
  public:
    /// @brief alias of default memory and host mirrorspace
    using DV = Kokkos::DualView<T*>;
    using MemorySpace = typename DV::t_dev::memory_space;
    using HostMirrorSpace = typename DV::t_host::memory_space;

    inline View(): __len(0) {}

    /// @brief Constructor
    inline void init(size_t len, size_t cap) {
      __len = len;
      this->realloc(cap);
    }
    inline void init(size_t len) { init(len, len); }

    /// @brief Modifiers
    inline void push_back(T in) { this->h_view(__len++) = in; }
    inline T pop_back() { return this->h_view(--__len); }

    /// @brief Capacity
    inline size_t size() const { return __len; }
    inline void resize(size_t n) { __len = n; }

    /// @brief Element access, for host
    inline T& operator[](int i) const { return this->h_view(i); }

    /// @brief Element access, for device
    KOKKOS_INLINE_FUNCTION
    T& operator()(int i) const { return this->d_view(i); }

    /// @brief The difference between 'erase' and 'remove' is 'remove' use value 
    ///     and 'erase' use pointer. 'erase' keeps the order while 'remove' don't.
    inline bool remove(T obj) {
      for (int i=0; i<__len; i++)
        if (this->h_view(i) == obj) {
          this->h_view(i) = this->h_view(--__len);
          return true;
        }
      return false;
    }

  private:
    size_t __len;
}; // class View


/**
 * @class Vector
 * @brief 3 dimensional vector
 * @details The following functions are operations of vectors.
 */
class Vector {
  public:
    KOKKOS_INLINE_FUNCTION Vector(): __data{0} {}
    KOKKOS_INLINE_FUNCTION
    Vector(double i, double j, double k): __data{i, j, k} {}
    KOKKOS_INLINE_FUNCTION
    double& operator[](int i) { return __data[i]; }
    KOKKOS_INLINE_FUNCTION
    double operator[](int i) const { return __data[i]; }

    /// @brief add, subtract and self-add
    KOKKOS_INLINE_FUNCTION
    Vector operator+(const Vector& p) const {
      return Vector(__data[0] + p[0], __data[1] + p[1], __data[2] + p[2]);
    }
    KOKKOS_INLINE_FUNCTION
    Vector operator-() const {
      return Vector(-__data[0], -__data[1], -__data[2]);
    }
    KOKKOS_INLINE_FUNCTION
    Vector operator-(const Vector& p) const {
      return Vector(__data[0] - p[0], __data[1] - p[1], __data[2] - p[2]);
    }
    KOKKOS_INLINE_FUNCTION
    Vector& operator+=(const Vector& p) {
      __data[0] += p[0]; __data[1] += p[1]; __data[2] += p[2];
      return *this;
    }

    /// @brief quantity product
    KOKKOS_INLINE_FUNCTION
    Vector operator*(double n) const {
      return Vector(n * __data[0], n * __data[1], n * __data[2]);
    }
    
    /// @brief 1/n times vector
    KOKKOS_INLINE_FUNCTION
    Vector operator/(double n) const {
      return Vector(__data[0] / n, __data[1] / n, __data[2] / n);
    }

    /// @brief dot product of two vectors
    KOKKOS_INLINE_FUNCTION
    double operator*(const Vector& p) const {
      return __data[0] * p[0] + __data[1] * p[1] + __data[2] * p[2];
    }

    /// @brief compare two vectors
    KOKKOS_INLINE_FUNCTION 
    bool operator==(const Vector& p) const {
      return __data[0]==p[0] && __data[1]==p[1] && __data[2]==p[2];
    }
    KOKKOS_INLINE_FUNCTION
    bool operator!=(const Vector& p) const {
      return __data[0]!=p[0] || __data[1]!=p[1] || __data[2]!=p[2];
    }

  private:
    double __data[3];
}; // class Vector

/// @brief quantity product for n on the left hand
KOKKOS_INLINE_FUNCTION
Vector operator*(double n, const Vector& p) { return p * n; }

/// @brief cross product of two vectors
KOKKOS_INLINE_FUNCTION
Vector cross(const Vector& p1, const Vector& p2) {
  return Vector( p1[1]*p2[2] - p1[2]*p2[1], p1[2]*p2[0] - p1[0]*p2[2],
    p1[0]*p2[1] - p1[1]*p2[0]);
}

/// @brief length of vector p
KOKKOS_INLINE_FUNCTION
double mod(const Vector& p) { return Kokkos::sqrt(p * p); }

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
    KOKKOS_INLINE_FUNCTION Pair(): __data{0} {}
    KOKKOS_INLINE_FUNCTION
    Pair(T a, T b): __data{a, b} {}
    KOKKOS_INLINE_FUNCTION
    T& operator[](int i) { return __data[i]; }
    KOKKOS_INLINE_FUNCTION
    T operator[](int i) const { return __data[i]; }
    
    /// @brief If node 'a' is connected by this bond, return the other node.
    ///     Return -1 else.
    KOKKOS_INLINE_FUNCTION
    T find(T a) const {
      return (__data[0] == a) ?
          __data[1] : ((__data[1] == a) ?
              __data[0] : -1);
    }

    /// @brief replace a with b，return false if failed，otherwise true.
    KOKKOS_INLINE_FUNCTION
    bool replace(T a, T b) {
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
    KOKKOS_INLINE_FUNCTION
    bool operator==(Pair<T> b) const {
      return (__data[0] == b[0] && __data[1] == b[1]) ||
          (__data[1] == b[0] && __data[0] == b[1]);
    }

    /// @brief add and substrct, for 2 dimensional coordinates
    KOKKOS_INLINE_FUNCTION
    Pair<T> operator-(Pair<T> b) const {
      return Pair<T>(__data[0] - b[0], __data[1] - b[1]);
    }
    KOKKOS_INLINE_FUNCTION
    Pair<T> operator+(Pair<T> b) const {
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
    KOKKOS_INLINE_FUNCTION
    Vector gen_vector(double end) const {
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

// reduction identity of Vector, must be defined in Kokkos namespace
namespace Kokkos {

template<>
struct reduction_identity<CoreMath::Vector> {
  KOKKOS_FORCEINLINE_FUNCTION
  static CoreMath::Vector sum() {
    return CoreMath::Vector();
  }
};

} // namespace Kokkos

#endif // QUADRATUBE_CORE_MATH_H_
