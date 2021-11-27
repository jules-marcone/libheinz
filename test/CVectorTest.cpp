#include "Vectors3D.h"
#include "Complex.h"
#include "Tests/GTestWrapper/google_test.h"

class CVectorTest : public ::testing::Test {
};

TEST_F(CVectorTest, TrivialOperations)
{
    R3 vec_k(1., 2., 3.);
    CHECK(complex_t(3., 0.) == vec_k.complex().z());

    C3 vec_c(complex_t(1., 3.), complex_t(2., -5.), complex_t(-3., -4.));
    CHECK(complex_t(3., 4.) == -vec_c.z());
    EXPECT_DOUBLE_EQ(vec_c.mag(), 8.);
}

TEST_F(CVectorTest, BasicArithmetics)
{
    // Dot product defined in Vec3
    C3 vec_a(complex_t(1., 0.), complex_t(2., 0.), complex_t(3., 0.));
    C3 vec_b(complex_t(2., 0.), complex_t(3., 0.), complex_t(4., 0.));
    C3 vec_c(complex_t(1., 1.), complex_t(2., -5.), complex_t(3., 4.));
    EXPECT_TRUE(vec_a.dot(vec_b) == complex_t(20., 0));
    EXPECT_TRUE(vec_a.dot(vec_c) == complex_t(14., 3.));
    EXPECT_TRUE(vec_c.dot(vec_b) == complex_t(20., -3.));
    EXPECT_TRUE(vec_a.dot(vec_a) == complex_t(14., 0));
    EXPECT_TRUE(vec_c.dot(vec_c) == complex_t(56., 0));

    // f = f_re + j*f_im
    C3 vec_e(1., 2., 3.);
    C3 vec_f(5., 6., 7.);
    CHECK(C3(complex_t(1., 5.), complex_t(2., 6), complex_t(3, 7)) == vec_e + I * vec_f);
}
