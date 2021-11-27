#include <heinz/Vectors3D.h>
#include "catch.hpp"

TEST_CASE("KVector:BasicMethods", "")
{
    R3 v;
    CHECK(v.x() == double(0));
    CHECK(v.y() == double(0));
    CHECK(v.z() == double(0));
    R3 v2(1., 2., 3.);
    CHECK(v2.x() == double(1.));
    CHECK(v2.y() == double(2.));
    CHECK(v2.z() == double(3.));
    v2.setX(10.);
    v2.setY(20.);
    v2.setZ(30.);
    CHECK(v2.x() == double(10.));
    CHECK(v2.y() == double(20.));
    CHECK(v2.z() == double(30.));
    v2 = {1., 2., 3.};
    CHECK(v2.x() == double(1.));
    CHECK(v2.y() == double(2.));
    CHECK(v2.z() == double(3.));

    R3 v3(1., 2., 3.);
    CHECK(v3.mag2() == 1 * 1 + 2 * 2 + 3 * 3);
    CHECK(v3.mag2() == v3.mag() * v3.mag());
    CHECK(v3.magxy2() == 1 * 1 + 2 * 2);
    CHECK(v3.magxy2() == Approx(v3.magxy() * v3.magxy()).epsilon(1e-15));
    CHECK(v3.magxy() == std::sqrt(1 * 1 + 2 * 2));
    CHECK(v3.mag() == std::sqrt(1 * 1 + 2 * 2 + 3 * 3));
}

TEST_CASE("KVector:BasicArithmetics", "")
{
    // assignment, self assignment, copy constructor
    R3 v1;
    R3 v2(v1);
    CHECK(v2.x() == double(0));
    CHECK(v2.y() == double(0));
    CHECK(v2.z() == double(0));
    v2 = {1., 2., 3.};
    CHECK(v2.x() == double(1));
    CHECK(v2.y() == double(2));
    CHECK(v2.z() == double(3));
    R3 v3(v2);
    CHECK(v3.x() == double(1));
    CHECK(v3.y() == double(2));
    CHECK(v3.z() == double(3));
    R3 v4 = v3;
    CHECK(v4.x() == double(1));
    CHECK(v4.y() == double(2));
    CHECK(v4.z() == double(3));
    // +=
    R3 a(1., 2., 3.);
    R3 b(10., 20., 30.);
    a += b;
    CHECK(a.x() == double(11));
    CHECK(a.y() == double(22));
    CHECK(a.z() == double(33));
    CHECK(b.x() == double(10));
    CHECK(b.y() == double(20));
    CHECK(b.z() == double(30));
    a = R3(1., 2., 3.);
    a += a;
    CHECK(a.x() == double(2.));
    CHECK(a.y() == double(4.));
    CHECK(a.z() == double(6.));
    // -=
    a = R3(1., 2., 3.);
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wself-assign-overloaded"
    // https://stackoverflow.com/a/57646029/1017348: a pretty stupid warning
    a -= a;
#pragma clang diagnostic pop
    CHECK(a.x() == double(0.));
    CHECK(a.y() == double(0.));
    CHECK(a.z() == double(0.));
    b = R3(1., 2., 3.);
    a -= b;
    CHECK(a.x() == double(-1.));
    CHECK(a.y() == double(-2.));
    CHECK(a.z() == double(-3.));
    // *= and /= with scalar
    a *= 10.0;
    CHECK(a.x() == double(-10.));
    CHECK(a.y() == double(-20.));
    CHECK(a.z() == double(-30.));
    a /= 10.;
    CHECK(a.x() == double(-1.));
    CHECK(a.y() == double(-2.));
    CHECK(a.z() == double(-3.));
    // unary minus
    a = R3(1., 2., 3.);
    b = -a;
    CHECK(a.x() == double(1.));
    CHECK(a.y() == double(2.));
    CHECK(a.z() == double(3.));
    CHECK(b.x() == double(-1.));
    CHECK(b.y() == double(-2.));
    CHECK(b.z() == double(-3.));
    // addition of two vector
    a = R3(1., 2., 3.);
    b = R3(10., 20., 30.);
    R3 c = a + b;
    CHECK(a.x() == double(1.));
    CHECK(a.y() == double(2.));
    CHECK(a.z() == double(3.));
    CHECK(b.x() == double(10.));
    CHECK(b.y() == double(20.));
    CHECK(b.z() == double(30.));
    CHECK(c.x() == double(11.));
    CHECK(c.y() == double(22.));
    CHECK(c.z() == double(33.));
    // substraction of two vectors
    c = b - a;
    CHECK(c.x() == double(9.));
    CHECK(c.y() == double(18.));
    CHECK(c.z() == double(27.));
    // multiplication by a scalar
    a = {1., 2., 3.};
    c = 2 * a * 2;
    CHECK(a.x() == double(1.));
    CHECK(a.y() == double(2.));
    CHECK(a.z() == double(3.));
    CHECK(c.x() == double(4.));
    CHECK(c.y() == double(8.));
    CHECK(c.z() == double(12.));
    // scalar product of two vectors
    a = {1., 2., 3.};
    b = {10., 10., 10.};
    CHECK(a.dot(b) == double(60));
    // crossproduct
    c = a.cross(b);
    CHECK(a.y() * b.z() - a.z() * b.y() == c.x());
    CHECK(a.z() * b.x() - a.x() * b.z() == c.y());
    CHECK(a.x() * b.y() - a.y() * b.x() == c.z());
    // equality
    a = {1., 2., 3.};
    CHECK(a == R3(1., 2., 3.));
    CHECK(a != R3(1., 1., 3.));
}
