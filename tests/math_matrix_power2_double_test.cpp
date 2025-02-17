#include <Eigen/Core>
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
    typedef std::chrono::high_resolution_clock Time;
    typedef std::chrono::milliseconds          ms;
    typedef std::chrono::duration<float>       fsec;
    auto                                       t0 = Time::now();
    callback();
    auto t1 = Time::now();
    fsec fs = t1 - t0;
    ms   d  = std::chrono::duration_cast<ms>(fs);
    return d.count();
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

template <class Type, std::size_t Size>
void
testMul()
{
    static_assert((Size & (Size - 1)) == 0, "Should be power of 2!");
    static const std::array<std::size_t, 10> arr   = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
    static const std::array<std::size_t, 10> iters = {1000000, 100000, 10000, 1000, 100, 10, 10, 1, 1, 1};

    using loc_type_strassen = matrix_strassen<Type, Size>;
    using loc_type_classic  = matrix_classic<Type, Size, Size>;
    using loc_type_eigen    = Eigen::Matrix<Type, Size, Size>;

    // Strassen
    auto              mAs = loc_type_strassen::get_rand_matrix();
    auto              mBs = loc_type_strassen::get_rand_matrix();
    loc_type_strassen mCs;
    // Classic
    auto             mAc = loc_type_classic::get_rand_matrix();
    auto             mBc = loc_type_classic::get_rand_matrix();
    loc_type_classic mCc;
    // Eigen
    auto mAe = loc_type_eigen::Random();
    auto mBe = loc_type_eigen::Random();

    int         i{};
    std::size_t iterations{};
    for (; i < arr.size(); i++) {
        if (Size == arr[i]) {
            break;
        }
    }
    if (i == arr.size()) {
        EXPECT_TRUE(false);
        return;
    } else {
        iterations = iters[i];
    }
    Type val{};

    auto strassen_time = measure([&]() {
        for (std::size_t k{}; k < iterations; k++) {
            mAs.mult(mBs, mCs);
            val = mCs.get(0, 0);
        }
    });
    std::cout << "Time of mult strassen " << Size << "x" << Size << " i=" << iterations << " checker: " << val << ": "
              << strassen_time << std::endl;

    auto classic_time = measure([&]() {
        for (std::size_t k{}; k < iterations; k++) {
            mAc.mult(mBc, mCc);
            val = mCc[0][0];
        }
    });
    std::cout << "Time of mult classic " << Size << "x" << Size << " i=" << iterations << " checker: " << val << ": "
              << classic_time << std::endl;

    auto eigen_time = measure([&]() {
        for (std::size_t k{}; k < iterations; k++) {
            val = (mAe * mBe)(0, 0);
        }
    });
    std::cout << "Time of mult eigen: " << Size << "x" << Size << " i=" << iterations << " checker: " << val << ": "
              << eigen_time << std::endl;

    auto strassen2_time = measure([&]() {
        for (std::size_t k{}; k < iterations; k++) {
            mAs.mult(mBs, mCs);
            mCs.add(mBs, mCs);
            val = mCs.get(0, 0);
        }
    });

    std::cout << "Time of mult + add strassen " << Size << "x" << Size << " i=" << iterations << " checker: " << val
              << ": " << strassen2_time << std::endl;

    auto classic2_time = measure([&]() {
        for (std::size_t k{}; k < iterations; k++) {
            mAc.mult(mBc, mCc);
            mCc.add(mBc, mCc);
            val = mCc[0][0];
        }
    });
    std::cout << "Time of mult + add classic " << Size << "x" << Size << " i=" << iterations << " checker: " << val
              << ": " << classic2_time << std::endl;

    auto eigen2_time = measure([&]() {
        for (std::size_t k{}; k < iterations; k++) {
            val = (mAe * mBe + mBe)(0, 0);
        }
    });
    std::cout << "Time of mult + add eigen: " << Size << "x" << Size << " i=" << iterations << " checker: " << val
              << ": " << eigen2_time << std::endl;
}

TEST(matrix_test, matrix_power_of_2_double)
{
    testMul<double, 2>();
    testMul<double, 4>();
    testMul<double, 8>();
    testMul<double, 16>();
    testMul<double, 32>();
    testMul<double, 64>();
    testMul<double, 128>();
    // testMul<double, 256>();
    // testMul<double, 512>();
    // testMul<double, 1024>();
}
