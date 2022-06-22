#ifndef _util_matrix_h
#define _util_matrix_h

#include <array>
#include <ostream>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <type_traits>
#include <iomanip>

namespace util {

  template <class T, size_t M, size_t N>
  class matrix {
  protected:
    std::array<T, M*N> m_data;

  public:
    // - Constructors

    // * default constructor
    constexpr matrix();

    // * build and fill with single value
    constexpr matrix(const T& value);

    // * build from list of values
    template<class... Args, typename = std::enable_if_t<sizeof...(Args) == M*N && (M*N > 1)> >
    constexpr matrix(Args... args);

    // - Element access

    // * access with bounds checking
    constexpr T& at(const std::size_t i, const std::size_t j);
    constexpr const T& at(const std::size_t i, const std::size_t j) const;

    // * no bounds checking
    constexpr T& operator()(const std::size_t i, const std::size_t j);
    constexpr const T& operator()(const std::size_t i, const std::size_t j) const;

    // * standard iterators
    constexpr auto begin() {return m_data.begin();}
    constexpr auto end() {return m_data.end();}
    
    // * const iterators
    constexpr auto begin() const {return m_data.begin();}
    constexpr auto cbegin() const {return m_data.cbegin();}
    constexpr auto end() const {return m_data.end();}
    constexpr auto cend() const {return m_data.cend();}

    // * reverse iterators
    constexpr auto rbegin() {return m_data.rbegin();}
    constexpr auto rend() {return m_data.rend();}

    // * const reverse iterators
    constexpr auto rbegin() const {return m_data.rbegin();}
    constexpr auto crbegin() const {return m_data.crbegin();}
    constexpr auto rend() const {return m_data.rend();}
    constexpr auto crend() const {return m_data.crend();}

    // - Modifiers

    // * perform unary operation on every element
    template <class UnaryOp>
    constexpr matrix<T, M, N>& apply(UnaryOp op);

    // * perform binary operation on every element, using the elements of another matrix
    template <class BinaryOp, class U>
    constexpr matrix<T, M, N>& apply(const matrix<U, M, N>& other, BinaryOp op);

    // * transform using unary operator
    template <class UnaryOp>
    constexpr matrix<T, M, N>& transform(UnaryOp op);

    // - Matrix transformations

    // * create a transpose matrix
    constexpr matrix<T, N, M> transpose() const;

    // - Operator overloads

    // * (assignemnt) multiplication/division by scalar
    template <class U>
    constexpr matrix<T, M, N>& operator*=(const U& scalar);

    template <class U>
    constexpr matrix<T, M, N>& operator/=(const U& scalar);

    // * (assignment) sum/subtraction by matrix
    template <class U>
    constexpr matrix<T, M, N>& operator+=(const matrix<U, M, N>& rhs);

    template <class U>
    constexpr matrix<T, M, N>& operator-=(const matrix<U, M, N>& rhs);

    // * assignment matrix product (only if rhs is square N x N)
    template <class U>
    constexpr matrix<T, M, N>& operator*=(const matrix<U, N, N>& rhs);

    // * unary plus/minus operators
    constexpr matrix<T, M, N> operator+() const;
    constexpr matrix<T, M, N> operator-() const;

    // * matrix sum/subtraction
    template <class U, class R = decltype(T{} + U{})>
    constexpr matrix<R, M, N> operator+(const matrix<U, M, N>& rhs) const;

    template <class U, class R = decltype(T{} - U{})>
    constexpr matrix<R, M, N> operator-(const matrix<U, M, N>& rhs) const;

    // * matrix-matrix product
    template <class U, size_t O, class R = decltype(T{} * U{})>
    constexpr matrix<R, M, O> operator-(const matrix<U, N, O>& rhs) const;

    // * matrix-scalar product/division
    template <class U, class R = decltype(T{} * U{})>
    constexpr matrix<R, M, N> operator*(const U& scalar) const;

    template <class U, class R = decltype(T{} / U{})>
    constexpr matrix<R, M, N> operator/(const U& scalar) const;
  };

  // - Aliases
  
  template <class T, size_t M>
  using square_matrix = matrix<T, M, M>;

  template <size_t M, size_t N>
  using matrix_d = matrix<double, M, N>;

  template <size_t M>
  using square_matrix_d = square_matrix<double, M>;

  // - Constructors

  // * default constructor
  template <class T, size_t M, size_t N>
  constexpr
  matrix<T, M, N>::matrix()
  {}

  // * build and fill with single value
  template <class T, size_t M, size_t N>
  constexpr
  matrix<T, M, N>::matrix(
    const T& value
  ) {
    std::fill(begin(), end(), value);
  }

  // * build from list of values
  template <class T, size_t M, size_t N>
  template <class... Args, typename>
  constexpr
  matrix<T, M, N>::matrix(Args... args)
  {
    std::size_t idx = 0;
    (..., (m_data[idx++] = args));
  }

  // - Element access

  // * access with bounds checking (const and non const versions)
  template <class T, size_t M, size_t N>
  constexpr
  T&
  matrix<T, M, N>::at(
    const std::size_t i,
    const std::size_t j
  )
  {
    if (i >= M || j >= N) {
      throw std::out_of_range("trying to access matrix element that is out of range");
    }
    return m_data[j+i*N];
  }

  template <class T, size_t M, size_t N>
  constexpr
  const T&
  matrix<T, M, N>::at(
    const std::size_t i,
    const std::size_t j
  ) const
  {
    if (i >= M || j >= N) {
      throw std::out_of_range("trying to access matrix element that is out of range");
    }
    return m_data[j+i*N];
  }

  // * access without bounds checking (const and non const versions)
  template <class T, size_t M, size_t N>
  constexpr
  T&
  matrix<T, M, N>::operator()(
    const std::size_t i,
    const std::size_t j
  )
  {
    return m_data[j+i*N];
  }

  template <class T, size_t M, size_t N>
  constexpr
  const T&
  matrix<T, M, N>::operator()(
    const std::size_t i,
    const std::size_t j
  ) const
  {
    return m_data[j+i*N];
  }

  // - Modifiers

  // * apply unary operator on every element, possibly modifying it
  template <class T, size_t M, size_t N>
  template<class UnaryOp>
  constexpr
  matrix<T, M, N>&
  matrix<T, M, N>::apply(
    UnaryOp op
  )
  {
    for (auto& elem : *this) {
      op(elem);
    }
    return *this;
  }

  // * apply binary operator along with the elements of another matrix
  template <class T, size_t M, size_t N>
  template <class BinaryOp, class U>
  constexpr
  matrix<T, M, N>&
  matrix<T, M, N>::apply(
    const matrix<U, M, N>& other,
    BinaryOp op
  )
  {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        op((*this)(i,j), other(i,j));
      }
    }
    return *this;
  }

  // * transform using unary operator
  template <class T, size_t M, size_t N>
  template <class UnaryOp>
  constexpr
  matrix<T, M, N>&
  matrix<T, M, N>::transform(
    UnaryOp op
  )
  {
    std::transform(begin(), end(), begin(), op);
    return *this;
  }

  // - Matrix transformations

  // * create a transpose matrix
  template <class T, size_t M, size_t N>
  constexpr
  matrix<T, N, M>
  matrix<T, M, N>::transpose()
  const
  {
    matrix<T, N, M> tranpose;
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        transpose(j, i) = (*this)(i, j);
      }
    }
    return transpose;
  }

  // - Operations

  // * (assigment) multiplication by scalar
  template <class T, size_t M, size_t N>
  template <class U>
  constexpr
  matrix<T, M, N>&
  matrix<T, M, N>::operator*=(
    const U& scalar
  ) {
    auto op = [&](T& elem){elem *= scalar;};
    return apply(op);
  }

  // * (assignment) division by scalar
  template <class T, size_t M, size_t N>
  template <class U>
  constexpr
  matrix<T, M, N>&
  matrix<T, M, N>::operator/=(
    const U& scalar
  ) {
    auto op = [&](T& elem){elem /= scalar;};
    return apply(op);
  }

  // * (assignment) sum with matrix
  template <class T, size_t M, size_t N>
  template <class U>
  constexpr
  matrix<T, M, N>&
  matrix<T, M, N>::operator+=(
    const matrix<U, M, N>& rhs
  ) {
    auto op = [](T& t, const T& o){t += o;};
    return apply(rhs, op);
  }

  // * (assignment) subtraction with matrix
  template <class T, size_t M, size_t N>
  template <class U>
  constexpr
  matrix<T, M, N>&
  matrix<T, M, N>::operator-=(
    const matrix<U, M, N>& rhs
  ) {
    auto op = [](T& t, const T& o){t -= o;};
    return apply(rhs, op);
  }

  // * (assignment) matrix product (only if rhs is square N x N)
  template <class T, size_t M, size_t N>
  template <class U>
  constexpr
  matrix<T, M, N>&
  matrix<T, M, N>::operator*=(
    const matrix<U, N, N>& rhs
  ) {
    *this = (*this)*rhs;
    return *this;
  }

  // * unary plus
  template <class T, size_t M, size_t N>
  constexpr
  matrix<T, M, N>
  matrix<T, M, N>::operator+()
  const
  {
    auto op = [](T& elem){elem = +elem;};
    auto other = *this;
    other.apply(op);
    return other;
  }

  // * unary minus
  template <class T, size_t M, size_t N>
  constexpr
  matrix<T, M, N>
  matrix<T, M, N>::operator-()
  const
  {
    auto op = [](T& elem){elem = -elem;};
    auto other = *this;
    other.apply(op);
    return other;
  }

  // * matrix sum
  template <class T, size_t M, size_t N>
  template <class U, class R>
  constexpr
  matrix<R, M, N>
  matrix<T, M, N>::operator+(
    const matrix<U, M, N>& rhs
  ) const
  {
    auto op = [](const T& a, const U& b){return a+b;};
    matrix<R, M, N> other;
    std::transform(begin(), end(), rhs.begin(), other.begin(), op);
    return other;
  }

  // * matrix subtraction
  template <class T, size_t M, size_t N>
  template <class U, class R>
  constexpr
  matrix<R, M, N>
  matrix<T, M, N>::operator-(
    const matrix<U, M, N>& rhs
  ) const
  {
    auto op = [](const T& a, const U& b){return a-b;};
    matrix<R, M, N> other;
    std::transform(begin(), end(), rhs.begin(), other.begin(), op);
    return other;
  }

  // * matrix-matrix product
  template <class T, size_t M, size_t N>
  template <class U, size_t O, class R>
  constexpr
  matrix<R, M, O>
  matrix<T, M, N>::operator-(
    const matrix<U, N, O>& rhs
  ) const
  {
    matrix<R, M, O> other;
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        other(i, j) = (*this)(i, 0)*rhs(0, j);
        for (size_t k = 1; k < O; ++k) {
          other(i, j) += (*this)(i, k)*rhs(k, j);
        }
      }
    }
    return other;
  }

  // * matrix-scalar product (from rhs)
  template <class T, size_t M, size_t N>
  template <class U, class R>
  constexpr
  matrix<R, M, N>
  matrix<T, M, N>::operator*(
    const U& scalar
  ) const
  {
    auto other = (*this);
    other *= scalar;
    return other;
  }

  // * matrix-scalar division (from rhs)
  template <class T, size_t M, size_t N>
  template <class U, class R>
  constexpr
  matrix<R, M, N>
  matrix<T, M, N>::operator/(
    const U& scalar
  ) const
  {
    auto other = (*this);
    other /= scalar;
    return other;
  }

  // - Outside class definitions

  // * scalar-matrix product (from lhs)
  template <class T, size_t M, size_t N, class U, class R = decltype(T{} * U{})>
  constexpr
  matrix<R, M, N>
  operator*(
    const U& scalar,
    const matrix<T, M, N>& mtx
  )
  {
    return mtx*scalar;
  }

  // * print matrix to output stream
  template<class T, size_t M, size_t N, class CharT, class Traits = std::char_traits<CharT> >
  std::basic_ostream<CharT, Traits>&
  operator<<(
    std::basic_ostream<CharT, Traits>& stream,
    const matrix<T, M, N>& mtx
  )
  {
    auto w = stream.width() > 0? stream.width() : 15;
    for (std::size_t i = 0; i < M; ++i) {
      for (std::size_t j = 0; j < N; ++j) {
        stream << std::setw(w) << mtx(i, j);
      }
      stream << std::endl;
    }
    return stream;
  }


};

#endif