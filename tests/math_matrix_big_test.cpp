#include <xitren/math/matrix.hpp>

#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace xitren::math;

TEST(matrix_test, matrix_func)
{
    using loc_type = matrix<double, 1024, 1025>;

    loc_type mA{};

    std::cout << "batch_value: " << loc_type::batch_value << std::endl;
    std::cout << "batch_rows: " << loc_type::batch_rows << std::endl;
    std::cout << "batch_columns: " << loc_type::batch_columns << std::endl;
    std::cout << "rest_rows: " << loc_type::rest_rows << std::endl;
    std::cout << "rest_columns: " << loc_type::rest_columns << std::endl;

    std::vector<double> init_vect{};
}
