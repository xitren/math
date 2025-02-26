#include <xitren/math/optimization.hpp>

#include <gtest/gtest.h>

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace xitren::math;

TEST(matrix_test, matrix_optimiation_1d)
{
    constexpr float precision{0.01};
    constexpr float result_val{4};

    optimization<float, 1, 1> simple(
        [](std::array<float, 1> in) -> std::array<float, 1> { return std::array<float, 1>{in[0] * in[0]}; },
        std::array<float, 1>{{16.}}, std::array<float, 1>{{2.}}, precision, 10s);

    simple.keep_running().wait(false);
    simple.thread().join();
    auto result = simple.result();
    for (std::size_t i{}; i < result.size(); i++) {
        std::cout << "value " << result[i] << " i " << i << std::endl;
    }
    EXPECT_TRUE(std::abs(result[0] - result_val) < precision);
}
