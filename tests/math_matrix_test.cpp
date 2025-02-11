#include <xitren/math/matrix.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace xitren::math;

template <class Type, std::size_t Size>
static void
print_matrix(matrix<Type, Size> const& pr)
{
    for (int i = 0; i < Size; i++) {
        for (int j = 0; j < Size; j++) {
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

    matrix<int, 2> mA{A};
    matrix<int, 2> mB{B};

    print_matrix(mA);
    print_matrix(mB);
    auto mC = mA * mB;
    print_matrix(mC);

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

    matrix<int, 4> mA{A};
    matrix<int, 4> mB{B};

    print_matrix(mA);
    print_matrix(mB);
    auto mC = mA * mB;
    print_matrix(mC);

    for (std::size_t l{}; l < 4; l++) {
        for (std::size_t m{}; m < 4; m++) {
            EXPECT_EQ(C[(l)*4 + m], mC.get(l, m));
        }
    }
}
