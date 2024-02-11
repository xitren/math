#include <xitren/math/bezier.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <iostream>

using bezier = xitren::math::bezier_cubic_d;

static bool
double_compare(const double val1, const double val2)
{
    static constexpr double precision = 0.0000001;
    if (val1 > val2) {
        return (val1 - val2 <= precision);
    } else {
        return (val2 - val1 <= precision);
    }
}

static bool
relative_compare(const double computed_val, const double real_val)
{
    static constexpr double relative_precision = 0.001;    // %

    if (real_val == 0.0f) {
        return false;
    }

    double error_ = std::fabs((computed_val - real_val) / real_val);

    if (error_ <= relative_precision) {
        return true;
    } else {
        return false;
    }
}

static bool
match(bezier::point_t p1, bezier::point_t p2)
{
    auto ret = double_compare(p1.x, p2.x) && double_compare(p1.y, p2.y);
    if (!ret) {
        std::cout << "p0(" << p1.x << ", " << p1.x << ") != p1(" << p2.y << ", " << p2.y << ")"
                  << std::endl;
    }
    return ret;
}

class bezier_cubic_manual {

public:
    bezier_cubic_manual() = default;

    bezier::point_t
    value(const bezier::points_array_t& set, const double t)
    {
        return {calc(set[0].x, set[1].x, set[2].x, set[3].x, t),
                calc(set[0].y, set[1].y, set[2].y, set[3].y, t)};
    }

private:
    bezier::coord_t
    calc(bezier::coord_t p0, bezier::coord_t p1, bezier::coord_t p2, bezier::coord_t p3,
         const double t)
    {
        p0 = p0 + (p1 - p0) * t;
        p1 = p1 + (p2 - p1) * t;
        p2 = p2 + (p3 - p2) * t;

        p0 = p0 + (p1 - p0) * t;
        p1 = p1 + (p2 - p1) * t;

        p0 = p0 + (p1 - p0) * t;

        return p0;
    }
};

TEST(bezier_test, base_point)
{
    bezier::point_t        p0{0, 0}, p1{0, 2}, p2{4, 2}, p3{4, 4};
    bezier::points_array_t base_points{p0, p1, p2, p3};

    bezier curve;
    curve.update(base_points);

    EXPECT_TRUE(match(curve.value(0), p0));
    EXPECT_TRUE(match(curve.value(1), p3));
}

TEST(bezier_test, alternative_calculate)
{
    bezier::point_t        p0{0, 0}, p1{0, 2}, p2{4, 2}, p3{4, 4};
    bezier::points_array_t base_points{p0, p1, p2, p3};

    bezier curve(base_points);

    bezier_cubic_manual curve2;

    double t{0};
    while (t < 1) {
        auto ret = match(curve.value(t), curve2.value(base_points, t));
        if (!ret) {
            std::cout << "Error: t = " << t << std::endl;
        }
        EXPECT_TRUE(ret);
        t += 0.005;
    }
}

TEST(bezier_test, curve_length)
{
    bezier::point_t        p0{0, 0}, p1{0, 2}, p2{4, 2}, p3{4, 4};
    bezier::points_array_t base_points{p0, p1, p2, p3};

    bezier curve(base_points);

    // compare with manual calculated value
    EXPECT_TRUE(relative_compare(curve.get_length(5), 5.9152));

    // compare with manual calculated value
    EXPECT_TRUE(relative_compare(curve.get_length(1), 5.6568));

    EXPECT_TRUE(double_compare(curve.get_length(0), 0));
}

TEST(bezier_test, curve_reveresed_length)
{
    bezier::point_t        p0{0, 0}, p1{0, -2}, p2{-4, -2}, p3{-4, -4};
    bezier::points_array_t base_points{p0, p1, p2, p3};

    bezier curve(base_points);

    // compare with manual calculated value
    EXPECT_TRUE(relative_compare(curve.get_length(5), 5.9152));
}

TEST(bezier_test, straight_y_length)
{
    bezier::point_t        p0{0, 0}, p1{0, 1}, p2{0, 2}, p3{0, 4};
    bezier::points_array_t base_points{p0, p1, p2, p3};

    bezier curve(base_points);

    EXPECT_TRUE(relative_compare(curve.get_length(5), 4));
}

TEST(bezier_test, straight_x_length)
{
    bezier::point_t        p0{0, 0}, p1{1, 0}, p2{2, 0}, p3{4, 0};
    bezier::points_array_t base_points{p0, p1, p2, p3};

    bezier curve(base_points);

    EXPECT_TRUE(relative_compare(curve.get_length(5), 4));
}
