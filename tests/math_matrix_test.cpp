#include <xitren/math/matrix_classic.hpp>
#include <xitren/math/matrix_strassen.hpp>

#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace xitren::math;

template <std::invocable<> Callable>
auto
measure(Callable callback)
{
    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::microseconds;
    auto start = high_resolution_clock::now();
    callback();
    auto end    = high_resolution_clock::now();
    auto ms_int = duration_cast<microseconds>(end - start);
    return ms_int;
}

template <class Type, std::size_t Size>
static void
print_matrix_strassen(matrix_strassen<Type, Size> const& pr)
{
    for (std::size_t i = 0; i < Size; i++) {
        for (std::size_t j = 0; j < Size; j++) {
            std::cout << pr.get(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

TEST(matrix_test, matrix_strassen_2x2)
{
    std::array<int, 4> A{{1, 3, 5, 7}};
    std::array<int, 4> B{{6, 8, 4, 2}};
    std::array<int, 4> C{{18, 14, 58, 54}};

    matrix_strassen<int, 2> mA{A};
    matrix_strassen<int, 2> mB{B};

    auto mC = mA * mB;

    EXPECT_EQ(C[0], mC.get(0, 0));
    EXPECT_EQ(C[1], mC.get(0, 1));
    EXPECT_EQ(C[2], mC.get(1, 0));
    EXPECT_EQ(C[3], mC.get(1, 1));
}

TEST(matrix_test, matrix_strassen_4x4)
{
    std::array<int, 16> A{{1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7}};
    std::array<int, 16> B{{7, 6, 5, 4, 3, 2, 1, 9, 8, 7, 6, 5, 4, 3, 2, 1}};
    std::array<int, 16> C{{53, 43, 33, 41, 141, 115, 89, 117, 94, 79, 64, 58, 119, 97, 75, 98}};

    matrix_strassen<int, 4> mA{A};
    matrix_strassen<int, 4> mB{B};

    auto mC = mA * mB;

    for (std::size_t l{}; l < 4; l++) {
        for (std::size_t m{}; m < 4; m++) {
            EXPECT_EQ(C[(l)*4 + m], mC.get(l, m));
        }
    }
}

constexpr std::size_t times_small = 1000000;
constexpr std::size_t times       = 10000;
constexpr std::size_t times_big   = 100;

TEST(matrix_test, matrix_strassen_2x2_mult_time)
{
    using loc_type = matrix_strassen<int, 4>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times_small; i++) {
        mCptr = std::make_shared<loc_type>(mA * mB);
    }
    std::cout << mCptr->get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_strassen_4x4_mult_time)
{
    using loc_type = matrix_strassen<int, 4>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times; i++) {
        mCptr = std::make_shared<loc_type>(mA * mB);
    }
    std::cout << mCptr->get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_strassen_8x8_mult_time)
{
    using loc_type = matrix_strassen<int, 8>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times; i++) {
        mCptr = std::make_shared<loc_type>(mA * mB);
    }
    std::cout << mCptr->get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_strassen_16x16_mult_time)
{
    using loc_type = matrix_strassen<int, 16>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times; i++) {
        mCptr = std::make_shared<loc_type>(mA * mB);
    }
    std::cout << mCptr->get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_strassen_32x32_mult_time)
{
    using loc_type = matrix_strassen<int, 32>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times_big; i++) {
        mCptr = std::make_shared<loc_type>(mA * mB);
    }
    std::cout << mCptr->get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_strassen_64x64_mult_time)
{
    using loc_type = matrix_strassen<int, 64>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times_big; i++) {
        mCptr = std::make_shared<loc_type>(mA * mB);
    }
    std::cout << mCptr->get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_strassen_128x128_mult_time)
{
    using loc_type = matrix_strassen<double, 128>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times_big; i++) {
        mCptr = std::make_shared<loc_type>(mA * mB);
    }
    std::cout << mCptr->get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_naive_2x2_mult_time)
{
    using loc_type = matrix<int, 2>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times_small; i++) {
        mCptr = std::make_shared<loc_type>(mA * mB);
    }
    std::cout << (*mCptr)[0][0] << std::endl;
}

TEST(matrix_test, matrix_naive_4x4_mult_time)
{
    using loc_type = matrix<int, 4>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times; i++) {
        mCptr = std::make_shared<loc_type>(mA * mB);
    }
    std::cout << (*mCptr)[0][0] << std::endl;
}

TEST(matrix_test, matrix_naive_8x8_mult_time)
{
    using loc_type = matrix<int, 8>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times; i++) {
        mCptr = std::make_shared<loc_type>(mA * mB);
    }
    std::cout << (*mCptr)[0][0] << std::endl;
}

TEST(matrix_test, matrix_naive_16x16_mult_time)
{
    using loc_type = matrix<int, 16>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times; i++) {
        mCptr = std::make_shared<loc_type>(mA * mB);
    }
    std::cout << (*mCptr)[0][0] << std::endl;
}

TEST(matrix_test, matrix_naive_32x32_mult_time)
{
    using loc_type = matrix<int, 32>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times_big; i++) {
        mCptr = std::make_shared<loc_type>(mA * mB);
    }
    std::cout << (*mCptr)[0][0] << std::endl;
}

TEST(matrix_test, matrix_naive_64x64_mult_time)
{
    using loc_type = matrix<int, 64>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times_big; i++) {
        mCptr = std::make_shared<loc_type>(mA * mB);
    }
    std::cout << (*mCptr)[0][0] << std::endl;
}

TEST(matrix_test, matrix_naive_128x128_mult_time)
{
    using loc_type = matrix<double, 128>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times_big; i++) {
        mCptr = std::make_shared<loc_type>(mA * mB);
    }
    std::cout << (*mCptr)[0][0] << std::endl;
}

TEST(matrix_test, matrix_strassen_2x2_add_time)
{
    using loc_type = matrix_strassen<int, 2>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times_small; i++) {
        mCptr = std::make_shared<loc_type>(mA + mB);
    }
    std::cout << mCptr->get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_strassen_4x4_add_time)
{
    using loc_type = matrix_strassen<int, 4>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times_small; i++) {
        mCptr = std::make_shared<loc_type>(mA + mB);
    }
    std::cout << mCptr->get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_strassen_8x8_add_time)
{
    using loc_type = matrix_strassen<int, 8>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times_small; i++) {
        mCptr = std::make_shared<loc_type>(mA + mB);
    }
    std::cout << mCptr->get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_strassen_16x16_add_time)
{
    using loc_type = matrix_strassen<int, 16>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times; i++) {
        mCptr = std::make_shared<loc_type>(mA + mB);
    }
    std::cout << mCptr->get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_strassen_32x32_add_time)
{
    using loc_type = matrix_strassen<int, 32>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times; i++) {
        mCptr = std::make_shared<loc_type>(mA + mB);
    }
    std::cout << mCptr->get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_strassen_64x64_add_time)
{
    using loc_type = matrix_strassen<int, 64>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times; i++) {
        mCptr = std::make_shared<loc_type>(mA + mB);
    }
    std::cout << mCptr->get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_strassen_128x128_add_time)
{
    using loc_type = matrix_strassen<double, 128>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times_big; i++) {
        mCptr = std::make_shared<loc_type>(mA + mB);
    }
    std::cout << mCptr->get(0, 0) << std::endl;
}

TEST(matrix_test, matrix_naive_2x2_add_time)
{
    using loc_type = matrix<int, 2>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times_small; i++) {
        mCptr = std::make_shared<loc_type>(mA + mB);
    }
    std::cout << (*mCptr)[0][0] << std::endl;
}

TEST(matrix_test, matrix_naive_4x4_add_time)
{
    using loc_type = matrix<int, 4>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times_small; i++) {
        mCptr = std::make_shared<loc_type>(mA + mB);
    }
    std::cout << (*mCptr)[0][0] << std::endl;
}

TEST(matrix_test, matrix_naive_8x8_add_time)
{
    using loc_type = matrix<int, 8>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times_small; i++) {
        mCptr = std::make_shared<loc_type>(mA + mB);
    }
    std::cout << (*mCptr)[0][0] << std::endl;
}

TEST(matrix_test, matrix_naive_16x16_add_time)
{
    using loc_type = matrix<int, 16>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times; i++) {
        mCptr = std::make_shared<loc_type>(mA + mB);
    }
    std::cout << (*mCptr)[0][0] << std::endl;
}

TEST(matrix_test, matrix_naive_32x32_add_time)
{
    using loc_type = matrix<int, 32>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times; i++) {
        mCptr = std::make_shared<loc_type>(mA + mB);
    }
    std::cout << (*mCptr)[0][0] << std::endl;
}

TEST(matrix_test, matrix_naive_64x64_add_time)
{
    using loc_type = matrix<int, 64>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times; i++) {
        mCptr = std::make_shared<loc_type>(mA + mB);
    }
    std::cout << (*mCptr)[0][0] << std::endl;
}

TEST(matrix_test, matrix_naive_128x128_add_time)
{
    using loc_type = matrix<double, 128>;

    auto mA = loc_type::get_rand_matrix();
    auto mB = loc_type::get_rand_matrix();

    auto mCptr = std::make_shared<loc_type>(mA * mB);
    for (std::size_t i{}; i < times_big; i++) {
        mCptr = std::make_shared<loc_type>(mA + mB);
    }
    std::cout << (*mCptr)[0][0] << std::endl;
}
