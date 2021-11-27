//  ************************************************************************************************
//
//  libheinz:  C++ base library of Heinz Maier-Leibnitz Zentrum
//
//! @file      Complex.h
//! @brief     Defines complex_t, and a few elementary functions
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libheinz
//! @license   Public Domain (BSD Zero Clause License, see LICENSE)
//! @copyright Forschungszentrum Jülich GmbH 2021
//! @authors   Scientific Computing Group at MLZ
//
//  ************************************************************************************************

#ifndef LIBHEINZ_COMPLEX_H
#define LIBHEINZ_COMPLEX_H

#include <complex>

using complex_t = std::complex<double>;

constexpr complex_t I = complex_t(0.0, 1.0);

//! Returns product I*z, where I is the imaginary unit.
inline complex_t mul_I(complex_t z)
{
    return complex_t(-z.imag(), z.real());
}

//! Returns exp(I*z), where I is the imaginary unit.
inline complex_t exp_I(complex_t z)
{
    return std::exp(complex_t(-z.imag(), z.real()));
}

#endif // LIBHEINZ_COMPLEX_H
