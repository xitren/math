#include <xitren/math/matrix.hpp>
#include <xitren/math/matrix_classic.hpp>

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

    static loc_type mA{{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}};
    static loc_type mB{{10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120}};

    print_matrix(mA);
    print_matrix(mB);

    static loc_type mC{};
    mA.add(mB, mC);

    print_matrix(mC);
    for (std::size_t l{}; l < 3; l++) {
        for (std::size_t m{}; m < 4; m++) {
            EXPECT_EQ(C[(l)*4 + m], mC.get(l, m));
        }
    }
}

TEST(matrix_big_test, matrix_func_mult_part)
{
    using loc_type_a = matrix<double, 2, 2, 2>;
    using loc_type_b = matrix<double, 2, 2, 2>;
    using loc_type_c = matrix<double, 2, 2, 2>;

    static loc_type_a mA1{{1, 2, 6, 7}};
    static loc_type_b mB1{{10, 20, 40, 50}};
    static loc_type_a mA2{{3, 4, 8, 9}};
    static loc_type_b mB2{{70, 80, 100, 110}};

    print_matrix(mA1);
    print_matrix(mB1);

    static loc_type_c mC1{};
    mA1.mult(mB1, mC1);
    print_matrix(mC1);

    static loc_type_c mC2{};
    mA2.mult(mB2, mC2);
    print_matrix(mC2);

    static loc_type_c mS{};
    mC1.mult(mC2, mS);
    print_matrix(mS);
}

TEST(matrix_big_test, matrix_func_mult)
{
    using loc_type_a = matrix<double, 3, 5, 2>;
    using loc_type_b = matrix<double, 5, 3, 2>;
    using loc_type_c = matrix<double, 3, 3, 2>;
    std::array<int, 9> C{{1350, 1500, 1650, 3100, 3500, 3900, 4850, 5500, 6150}};

    static loc_type_a mA{{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}};
    static loc_type_b mB{{10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150}};

    print_matrix(mA);
    print_matrix(mB);

    static loc_type_c mC{};
    mA.mult(mB, mC);

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

constexpr std::size_t times_big   = 100;
constexpr std::size_t times_small = 10000;

TEST(matrix_test, matrix_hybrid_64x64_mult_time)
{
    using loc_type1    = matrix<double, 64, 128, 64>;
    using loc_type2    = matrix<double, 128, 64, 64>;
    using loc_type_ret = matrix<double, 64, 64, 64>;

    static loc_type1 mA = loc_type1::get_rand_matrix();
    static loc_type2 mB = loc_type2::get_rand_matrix();

    static loc_type_ret mC{};
    mA.mult(mB, mC);
    for (std::size_t i{}; i < times_big; i++) {
        mA.mult(mB, mC);
    }
    std::cout << mC.get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_naive_64x64_mult_time)
{
    using loc_type1    = matrix_classic<double, 64, 128>;
    using loc_type2    = matrix_classic<double, 128, 64>;
    using loc_type_ret = matrix_classic<double, 64, 64>;

    static auto mA = loc_type1::get_rand_matrix();
    static auto mB = loc_type2::get_rand_matrix();

    static loc_type_ret mC{};
    mA.mult(mB, mC);
    for (std::size_t i{}; i < times_big; i++) {
        mA.mult(mB, mC);
    }
    std::cout << mC[0][0] << std::endl;
}

TEST(matrix_test, matrix_hybrid_77x77_mult_time)
{
    using loc_type1    = matrix<double, 77, 30, 16>;
    using loc_type2    = matrix<double, 30, 77, 16>;
    using loc_type_ret = matrix<double, 77, 77, 16>;

    static auto mA = loc_type1::get_rand_matrix();
    static auto mB = loc_type2::get_rand_matrix();

    static loc_type_ret mC{};
    mA.mult(mB, mC);
    for (std::size_t i{}; i < times_big; i++) {
        mA.mult(mB, mC);
    }
    std::cout << mC.get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_naive_77x77_mult_time)
{
    using loc_type1    = matrix_classic<double, 77, 30>;
    using loc_type2    = matrix_classic<double, 30, 77>;
    using loc_type_ret = matrix_classic<double, 77, 77>;

    static auto mA = loc_type1::get_rand_matrix();
    static auto mB = loc_type2::get_rand_matrix();

    static loc_type_ret mC{};
    mA.mult(mB, mC);
    for (std::size_t i{}; i < times_big; i++) {
        mA.mult(mB, mC);
    }
    std::cout << mC[0][0] << std::endl;
}

TEST(matrix_test, matrix_hybrid_64x64_mult_unit8_time)
{
    using loc_type1    = matrix<std::uint8_t, 64, 128, 64>;
    using loc_type2    = matrix<std::uint8_t, 128, 64, 64>;
    using loc_type_ret = matrix<std::uint8_t, 64, 64, 64>;

    static auto mA = loc_type1::get_rand_matrix();
    static auto mB = loc_type2::get_rand_matrix();

    static loc_type_ret mC{};
    mA.mult(mB, mC);
    for (std::size_t i{}; i < times_big; i++) {
        mA.mult(mB, mC);
    }
    std::cout << mC.get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_naive_64x64_mult_unit8_time)
{
    using loc_type1    = matrix_classic<std::uint8_t, 64, 128>;
    using loc_type2    = matrix_classic<std::uint8_t, 128, 64>;
    using loc_type_ret = matrix_classic<uint8_t, 64, 64>;

    static auto mA = loc_type1::get_rand_matrix();
    static auto mB = loc_type2::get_rand_matrix();

    static loc_type_ret mC{};
    mA.mult(mB, mC);
    for (std::size_t i{}; i < times_big; i++) {
        mA.mult(mB, mC);
    }
    std::cout << mC[0][0] << std::endl;
}

TEST(matrix_test, matrix_hybrid_77x77_mult_unit8_time)
{
    using loc_type1    = matrix<std::uint8_t, 77, 30, 16>;
    using loc_type2    = matrix<std::uint8_t, 30, 77, 16>;
    using loc_type_ret = matrix<std::uint8_t, 77, 77, 16>;

    static auto mA = loc_type1::get_rand_matrix();
    static auto mB = loc_type2::get_rand_matrix();

    static loc_type_ret mC{};
    mA.mult(mB, mC);
    for (std::size_t i{}; i < times_big; i++) {
        mA.mult(mB, mC);
    }
    std::cout << mC.get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_naive_77x77_mult_unit8_time)
{
    using loc_type1    = matrix_classic<std::uint8_t, 77, 30>;
    using loc_type2    = matrix_classic<std::uint8_t, 30, 77>;
    using loc_type_ret = matrix_classic<uint8_t, 77, 77>;

    static auto mA = loc_type1::get_rand_matrix();
    static auto mB = loc_type2::get_rand_matrix();

    static loc_type_ret mC{};
    mA.mult(mB, mC);
    for (std::size_t i{}; i < times_big; i++) {
        mA.mult(mB, mC);
    }
    std::cout << mC[0][0] << std::endl;
}

TEST(matrix_test, matrix_hybrid_64x128_add_time)
{
    using loc_type1    = matrix<double, 64, 128, 64>;
    using loc_type2    = matrix<double, 64, 128, 64>;
    using loc_type_ret = matrix<double, 64, 128, 64>;

    static auto mA = loc_type1::get_rand_matrix();
    static auto mB = loc_type2::get_rand_matrix();

    static loc_type_ret mC{};
    mA.add(mB, mC);
    for (std::size_t i{}; i < times_small; i++) {
        mA.add(mB, mC);
    }
    std::cout << mC.get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_naive_64x128_add_time)
{
    using loc_type1    = matrix_classic<double, 64, 128>;
    using loc_type2    = matrix_classic<double, 64, 128>;
    using loc_type_ret = matrix_classic<double, 64, 128>;

    static auto mA = loc_type1::get_rand_matrix();
    static auto mB = loc_type2::get_rand_matrix();

    static loc_type_ret mC{};
    mA.add(mB, mC);
    for (std::size_t i{}; i < times_small; i++) {
        mA.add(mB, mC);
    }
    std::cout << mC[0][0] << std::endl;
}

TEST(matrix_test, matrix_hybrid_77x77_add_time)
{
    using loc_type1    = matrix<double, 77, 77, 16>;
    using loc_type2    = matrix<double, 77, 77, 16>;
    using loc_type_ret = matrix<double, 77, 77, 16>;

    static auto mA = loc_type1::get_rand_matrix();
    static auto mB = loc_type2::get_rand_matrix();

    static loc_type_ret mC{};
    mA.add(mB, mC);
    for (std::size_t i{}; i < times_small; i++) {
        mA.add(mB, mC);
    }
    std::cout << mC.get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_naive_77x77_add_time)
{
    using loc_type1    = matrix_classic<double, 77, 77>;
    using loc_type2    = matrix_classic<double, 77, 77>;
    using loc_type_ret = matrix_classic<double, 77, 77>;

    static auto mA = loc_type1::get_rand_matrix();
    static auto mB = loc_type2::get_rand_matrix();

    static loc_type_ret mC{};
    mA.add(mB, mC);
    for (std::size_t i{}; i < times_small; i++) {
        mA.add(mB, mC);
    }
    std::cout << mC[0][0] << std::endl;
}
