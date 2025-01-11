#include <xitren/math/bezier.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <iostream>

using namespace xitren::math;

TEST(bezier_test, base_point)
{
    bezier_point<int>                p0{0, 0}, p1{0, 2}, p2{4, 2}, p3{4, 4};
    std::array<bezier_point<int>, 4> base_points{p0, p1, p2, p3};

    bezier_quadratic<int, 100> curve{};
    curve.update(base_points);

    EXPECT_TRUE(curve[0] == p0);
    EXPECT_TRUE(curve[99] == p3);
}

TEST(bezier_test, base_point2)
{
    bezier_point<int>                p0{-7, 7}, p1{-7, 7}, p2{7, 7}, p3{7, 7};
    std::array<bezier_point<int>, 4> base_points{p0, p1, p2, p3};

    bezier_quadratic<int, 100> curve{};
    curve.update(base_points);

    // EXPECT_TRUE(curve[0] == p0);
    EXPECT_TRUE(curve[99] == p3);
}
