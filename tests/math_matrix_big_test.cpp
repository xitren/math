#include <xitren/math/matrix.hpp>

#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace xitren::math;

template <class Type, std::size_t Rows, std::size_t Columns, std::size_t Batch>
static void
print_matrix(matrix<Type, Rows, Columns, Batch>& pr)
{
    for (std::size_t i = 0; i < Rows; i++) {
        for (std::size_t j = 0; j < Columns; j++) {
            auto a = pr.get(i, j);
            std::cout << a << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

TEST(matrix_big_test, matrix_func_add)
{
    using loc_type = matrix<double, 3, 4>;
    std::array<int, 12> C{{11, 22, 33, 44, 55, 66, 77, 88, 99, 110, 121, 132}};

    loc_type mA{{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}};
    loc_type mB{{10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120}};

    print_matrix(mA);
    print_matrix(mB);

    auto mC = mA + mB;

    print_matrix(mC);
    for (std::size_t l{}; l < 3; l++) {
        for (std::size_t m{}; m < 4; m++) {
            EXPECT_EQ(C[(l)*4 + m], mC.get(l, m));
        }
    }
}

TEST(matrix_big_test, matrix_func_mult_part)
{
    using loc_type_a = matrix<double, 2, 2>;
    using loc_type_b = matrix<double, 2, 2>;

    loc_type_a mA1{{1, 2, 6, 7}};
    loc_type_b mB1{{10, 20, 40, 50}};
    loc_type_a mA2{{3, 4, 8, 9}};
    loc_type_b mB2{{70, 80, 100, 110}};

    print_matrix(mA1);
    print_matrix(mB1);

    auto mC1 = mA1 * mB1;
    print_matrix(mC1);

    auto mC2 = mA2 * mB2;
    print_matrix(mC2);

    auto mS = mC1 + mC2;
    print_matrix(mS);
}

TEST(matrix_big_test, matrix_func_mult)
{
    using loc_type_a = matrix<double, 3, 5>;
    using loc_type_b = matrix<double, 5, 3>;
    std::array<int, 9> C{{1350, 1500, 1650, 3100, 3500, 3900, 4850, 5500, 6150}};

    loc_type_a mA{{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}};
    loc_type_b mB{{10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150}};

    print_matrix(mA);
    print_matrix(mB);

    auto mC = mA * mB;

    print_matrix(mC);
    for (std::size_t l{}; l < 3; l++) {
        for (std::size_t m{}; m < 3; m++) {
            EXPECT_EQ(C[(l)*3 + m], mC.get(l, m));
        }
    }
}

TEST(matrix_big_test, matrix_func2)
{
    using loc_type = matrix<double, 1025, 1025>;

    static loc_type mA{};

    std::cout << "batch_value: " << loc_type::batch_value << std::endl;
    std::cout << "batch_rows: " << loc_type::batch_rows << std::endl;
    std::cout << "batch_columns: " << loc_type::batch_columns << std::endl;
    std::cout << "rest_rows: " << loc_type::rest_rows << std::endl;
    std::cout << "rest_columns: " << loc_type::rest_columns << std::endl;

    EXPECT_EQ(loc_type::batch_value, 128);
    EXPECT_EQ(loc_type::batch_rows, 8);
    EXPECT_EQ(loc_type::batch_columns, 8);
    EXPECT_EQ(loc_type::rest_rows, 1);
    EXPECT_EQ(loc_type::rest_columns, 1);
}
