//  ************************************************************************************************
//
//  libheinz:  C++ base library of Heinz Maier-Leibnitz Zentrum
//
//! @file      heinz/Vectors3D.h
//! @brief     Defines and implements three-dimensional vector types I3, R3, C3
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libheinz
//! @license   Public Domain (BSD Zero Clause License, see LICENSE)
//! @copyright Forschungszentrum Jülich GmbH 2021
//! @authors   Scientific Computing Group at MLZ
//
//  ************************************************************************************************

#ifndef LIBHEINZ_VECTORS3D_H
#define LIBHEINZ_VECTORS3D_H

#include <complex>

template <class T> class Vec3;
using I3 = Vec3<int>;
using R3 = Vec3<double>;
using C3 = Vec3<std::complex<double>>;


//! Three-dimensional vector class, templated on base type.

//! Base type T can be int, double, or complex<double>.
//! The resulting classes shall be abbreviated I3, R3, C3.

template <class T> class Vec3 {
private:
    T m_x;
    T m_y;
    T m_z;

public:
    // -------------------------------------------------------------------------
    // Constructors and other set functions
    // -------------------------------------------------------------------------

    //! Constructs a vector from cartesian components.
    Vec3(const T x_, const T y_, const T z_) : m_x(x_), m_y(y_), m_z(z_) {}

    //! Constructs the null vector.
    Vec3() : Vec3{0, 0, 0} {}


    // -------------------------------------------------------------------------
    // Component access
    // -------------------------------------------------------------------------

    //! Returns x-component in cartesian coordinate system.
    inline T x() const { return m_x; }
    //! Returns y-component in cartesian coordinate system.
    inline T y() const { return m_y; }
    //! Returns z-component in cartesian coordinate system.
    inline T z() const { return m_z; }

    //! Sets x-component in cartesian coordinate system.
    void setX(const T& a) { m_x = a; }
    //! Sets y-component in cartesian coordinate system.
    void setY(const T& a) { m_y = a; }
    //! Sets z-component in cartesian coordinate system.
    void setZ(const T& a) { m_z = a; }

    // -------------------------------------------------------------------------
    // In-place operations
    // -------------------------------------------------------------------------

    //! Adds other vector to this, and returns result.
    Vec3<T>& operator+=(const Vec3<T>& v)
    {
        m_x += v.x();
        m_y += v.y();
        m_z += v.z();
        return *this;
    }

    //! Subtracts other vector from this, and returns result.
    Vec3<T>& operator-=(const Vec3<T>& v)
    {
        m_x -= v.x();
        m_y -= v.y();
        m_z -= v.z();
        return *this;
    }

#ifndef SWIG
    //! Multiplies this with a scalar, and returns result.
    template <class U> auto operator*=(U a)
    {
        m_x *= a;
        m_y *= a;
        m_z *= a;
        return *this;
    }
#endif // SWIG

#ifndef SWIG
    //! Divides this by a scalar, and returns result.
    template <class U> auto operator/=(U a)
    {
        m_x /= a;
        m_y /= a;
        m_z /= a;
        return *this;
    }
#endif // SWIG

    // -------------------------------------------------------------------------
    // Functions of this (with no further argument)
    // -------------------------------------------------------------------------

    //! Returns complex conjugate vector
    Vec3<T> conj() const;

    //! Returns magnitude squared of the vector.
    double mag2() const { return std::norm(x()) + std::norm(y()) + std::norm(z()); }

    //! Returns magnitude of the vector.
    double mag() const { return sqrt(mag2()); }

    //! Returns squared distance from z axis.
    double magxy2() const { return std::norm(x()) + std::norm(y()); }

    //! Returns distance from z axis.
    double magxy() const { return sqrt(magxy2()); }

    //! Returns unit vector in direction of this. Throws for null vector.
    Vec3<T> unit() const;

    //! Returns this, trivially converted to complex type.
    C3 complex() const;

    //! Returns real parts.
    R3 real() const;

    // -------------------------------------------------------------------------
    // Functions of this and another vector
    // -------------------------------------------------------------------------

    inline bool operator==(const Vec3<T>& other) const {
        return x()==other.x() && y()==other.y() && z()==other.z(); }

    inline bool operator!=(const Vec3<T>& other) const { return !(*this==other); }

#ifndef SWIG
    //! Returns dot product of vectors (antilinear in the first [=self] argument).
    template <class U> auto dot(const Vec3<U>& v) const;
#endif // SWIG

#ifndef SWIG
    //! Returns cross product of vectors (linear in both arguments).
    template <class U> auto cross(const Vec3<U>& v) const;
#endif // SWIG

    //! Returns projection of this onto other vector: (this*v)*v/|v|^2.
    Vec3<T> project(const Vec3<T>& v) const { return dot(v) * v / v.mag2(); }

    // -------------------------------------------------------------------------
    // Rotations
    // -------------------------------------------------------------------------

    // //! Returns result of rotation around x-axis.
    // Vec3<T> rotatedX(double a) const;
    //! Returns result of rotation around y-axis.
    Vec3<T> rotatedY(double a) const
    {
        return Vec3<T>(cos(a) * x() + sin(a) * z(), y(), -sin(a) * x() + cos(a) * z());
    }
    //! Returns result of rotation around z-axis.
    Vec3<T> rotatedZ(double a) const
    {
        return Vec3<T>(cos(a) * x() + sin(a) * y(), -sin(a) * x() + cos(a) * y(), z());
    }
    // //! Returns result of rotation around the axis specified by another vector.
    // Vec3<T> rotated(double a, const Vec3<T>& v) const;
};

#ifndef USER_API

// =============================================================================
// Non-member functions
// =============================================================================

//! Output to stream.
//! @relates Vec3
template <class T> std::ostream& operator<<(std::ostream& os, const Vec3<T>& a)
{
    return os << "(" << a.x() << "," << a.y() << "," << a.z() << ")";
}

// -----------------------------------------------------------------------------
// Unary operators
// -----------------------------------------------------------------------------

//! Unary plus.
//! @relates Vec3
template <class T> inline Vec3<T> operator+(const Vec3<T>& v)
{
    return v;
}

//! Unary minus.
//! @relates Vec3
template <class T> inline Vec3<T> operator-(const Vec3<T>& v)
{
    return {-v.x(), -v.y(), -v.z()};
}

// -----------------------------------------------------------------------------
// Binary operators
// -----------------------------------------------------------------------------

//! Addition of two vectors.
//! @relates Vec3
template <class T> inline Vec3<T> operator+(const Vec3<T>& a, const Vec3<T>& b)
{
    return {a.x() + b.x(), a.y() + b.y(), a.z() + b.z()};
}

//! Subtraction of two vectors.
//! @relates Vec3
template <class T> inline Vec3<T> operator-(const Vec3<T>& a, const Vec3<T>& b)
{
    return {a.x() - b.x(), a.y() - b.y(), a.z() - b.z()};
}

#ifndef SWIG
//! Multiplication vector by scalar.
//! @relates Vec3
template <class T, class U> inline auto operator*(const Vec3<T>& v, const U a)
{
    return Vec3<decltype(v.x() * a)>{v.x() * a, v.y() * a, v.z() * a};
}
#endif // SWIG

#ifndef SWIG
//! Multiplication scalar by vector.
//! @relates Vec3
template <class T, class U> inline auto operator*(const U a, const Vec3<T>& v)
{
    return Vec3<decltype(a * v.x())>{a * v.x(), a * v.y(), a * v.z()};
}
#endif // SWIG

// vector*vector not supported
//    (We do not provide the operator form a*b of the dot product:
//     Though nice to write, and in some cases perfectly justified,
//     in general it tends to make expressions more difficult to read.)

//! Division vector by scalar.
//! @relates Vec3
template <class T, class U> inline Vec3<T> operator/(const Vec3<T>& v, U a)
{
    return {v.x() / a, v.y() / a, v.z() / a};
}

// =============================================================================
// ?? for API generation ??
// =============================================================================

#ifndef SWIG
//! Returns dot product of (complex) vectors (antilinear in the first [=self] argument).
template <class T> template <class U> inline auto Vec3<T>::dot(const Vec3<U>& v) const
{
    Vec3<T> left_star = this->conj();
    return left_star.x() * v.x() + left_star.y() * v.y() + left_star.z() * v.z();
}
#endif // SWIG

#ifndef SWIG
//! Returns cross product of (complex) vectors.
template <class T> template <class U> inline auto Vec3<T>::cross(const Vec3<U>& v) const
{
    return Vec3<decltype(this->x() * v.x())>(y() * v.z() - v.y() * z(), z() * v.x() - v.z() * x(),
                                             x() * v.y() - v.x() * y());
}
#endif // SWIG

// -----------------------------------------------------------------------------
// Functions of this (with no further argument)
// -----------------------------------------------------------------------------

template <> inline R3 R3::conj() const
{
    return *this;
}

template <> inline C3 C3::conj() const
{
    return {std::conj(x()), std::conj(y()), std::conj(z())};
}

template <> inline C3 R3::complex() const
{
    return {x(), y(), z()};
}

template <> inline C3 C3::complex() const
{
    return *this;
}

template <> inline R3 R3::real() const
{
    return *this;
}

template <> inline R3 C3::real() const
{
    return {x().real(), y().real(), z().real()};
}

template <> inline R3 R3::unit() const
{
    double len = mag();
    if (len == 0.0)
        throw std::runtime_error("Cannot normalize zero vector");
    return {x() / len, y() / len, z() / len};
}

template <> inline C3 C3::unit() const
{
    double len = mag();
    if (len == 0.0)
        throw std::runtime_error("Cannot normalize zero vector");
    return {x() / len, y() / len, z() / len};
}

// -----------------------------------------------------------------------------
// Combine two vectors
// -----------------------------------------------------------------------------

namespace R3Util {

//! Returns polar angle.
inline double theta(const R3& a)
{
    return a.x() == 0.0 && a.y() == 0.0 && a.z() == 0.0 ? 0.0 : std::atan2(a.magxy(), a.z());
}

inline double phi(const R3& a)
{
    return a.x() == 0.0 && a.y() == 0.0 ? 0.0 : std::atan2(-a.y(), a.x());
}

inline double cosTheta(const R3& a)
{
    return a.mag() == 0 ? 1 : a.z() / a.mag();
}

inline double sin2Theta(const R3& a)
{
    return a.mag2() == 0 ? 0 : a.magxy2() / a.mag2();
}

inline double angle(const R3& a, const R3& b)
{
    double cosa = 0;
    double ptot = a.mag() * b.mag();
    if (ptot > 0) {
        cosa = a.dot(b) / ptot;
        if (cosa > 1)
            cosa = 1;
        if (cosa < -1)
            cosa = -1;
    }
    return std::acos(cosa);
}

} // namespace R3Util

#endif // USER_API
#endif // LIBHEINZ_VECTORS3D_H
