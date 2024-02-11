#include <xitren/math/kht_opt.hpp>

#include <gtest/gtest.h>

#include <iostream>

using namespace xitren::math;

template <size_t Size>
bool
arrays_match(const std::array<std::uint8_t, Size>& expected,
             const std::array<std::uint8_t, Size>& actual)
{
    for (size_t i{0}; i < Size; ++i) {
        if (expected[i] != actual[i]) {
            std::cout << "array[" << i << "] (" << static_cast<int>(actual[i]) << ") != expected["
                      << i << "] (" << static_cast<int>(expected[i]) << ")";
            return false;
        }
    }
    return true;
}

TEST(kht_test, base_test_opt_2)
{
    constexpr std::size_t width  = 2;
    constexpr std::size_t height = 2;
    vault<width, height>  image{{0, 1, 1, 0}};
    kht<width, height, std::uint8_t, true>::convert(image);
    std::array<std::uint8_t, width * height> result{1, 1, 0, 1};
    std::array<std::uint8_t, width * height> result_mirror{1, 0, 2, 0};
    EXPECT_TRUE(arrays_match(result, image.image()));
    EXPECT_TRUE(arrays_match(result_mirror, image.mirror()));
}

TEST(kht_test, base_test_opt_4)
{
    constexpr std::size_t width  = 4;
    constexpr std::size_t height = 4;
    vault<width, height>  image{{0, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 1, 0, 0}};
    kht<width, height, std::uint8_t, true>::convert(image);
    std::array<std::uint8_t, width * height> result{0, 2, 3, 0, 2, 0, 3, 0, 1, 2, 1, 0, 0, 2, 1, 0};
    std::array<std::uint8_t, width * height> result_mirror{0, 2, 3, 0, 0, 5, 0, 0,
                                                           2, 3, 0, 0, 3, 1, 0, 0};
    image.print();
    EXPECT_TRUE(arrays_match(result, image.image()));
    EXPECT_TRUE(arrays_match(result_mirror, image.mirror()));
}

std::array<std::array<std::uint8_t, 8 * 8>, 8> angles_test{
    {{
         1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
         0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
     },
     {
         1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
         0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
     },
     {
         1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
         0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
     },
     {
         1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
         0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0,
         0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
     },
     {
         1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
         0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
     },
     {
         1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
         1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
     },
     {
         1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
         1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
     },
     {
         1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
         0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
     }}};

std::array<std::array<std::uint8_t, 8 * 8>, 8> angles_test_mirror{
    {{
         1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
         0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
     },
     {
         0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
         0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
     },
     {
         0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
         0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
     },
     {
         0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
     },
     {
         0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
     },
     {
         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0,
         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
     },
     {
         0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
     },
     {
         0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
         0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
     }}};

TEST(kht_test, base_test_opt_8_batch_16)
{
    constexpr std::size_t                    width  = 8;
    constexpr std::size_t                    height = 8;
    std::array<std::uint8_t, width * height> test_zero{{}};
    auto                                     i{0};
    auto                                     i_mirror{0};
    for (auto test_i : angles_test) {
        vault<width, height> image{test_i};
        kht<width, height, std::uint16_t>::convert(image);
        auto i_max
            = std::max_element(image.image().begin(), image.image().end()) - image.image().begin();
        auto x = i_max % width;
        auto y = i_max / width;
        std::cout << "x = " << x << "; y = " << y << ";" << std::endl;
        EXPECT_TRUE((x == 0) && (y == i++));
    }
    for (auto test_i : angles_test_mirror) {
        vault<width, height> image{test_i};
        kht<width, height, std::uint16_t>::convert(image);
        auto i_max = std::max_element(image.mirror().begin(), image.mirror().end())
                     - image.mirror().begin();
        auto x = i_max % width;
        auto y = i_max / width;
        std::cout << "Mirror x = " << x << "; y = " << y << ";" << std::endl;
        EXPECT_TRUE((x == 0) && (y == i_mirror++));
    }
}

TEST(kht_test, base_test_opt_16_vertical_bold)
{
    constexpr std::size_t width  = 16;
    constexpr std::size_t height = 16;
    vault<width, height>  image{
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
         1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0}};
    kht<width, height>::convert(image);
    std::array<std::uint8_t, width * height> result{
        0,  0, 0, 0, 0,  0,  0, 0, 4, 8,  12, 12, 0, 0, 0, 0,  0,  0, 0, 0, 0,  0,  0, 0, 4, 12,
        16, 5, 0, 0, 0,  0,  0, 0, 0, 0,  0,  0,  0, 0, 8, 16, 13, 1, 0, 0, 0,  0,  0, 0, 0, 0,
        0,  0, 0, 0, 16, 13, 9, 0, 0, 0,  0,  0,  0, 0, 0, 0,  0,  0, 0, 8, 15, 12, 5, 0, 0, 0,
        0,  0, 0, 0, 0,  0,  0, 0, 4, 11, 12, 9,  4, 0, 0, 0,  0,  0, 0, 0, 0,  0,  0, 2, 7, 11,
        12, 7, 2, 0, 0,  0,  0, 0, 0, 0,  0,  0,  2, 5, 7, 10, 9,  6, 2, 0, 0,  0,  0, 0, 0, 0,
        0,  1, 4, 7, 9,  11, 7, 4, 1, 0,  0,  0,  0, 0, 0, 0,  1,  4, 6, 6, 8,  8,  7, 4, 1, 0,
        0,  0, 0, 0, 0,  1,  3, 4, 5, 8,  8,  7,  5, 3, 1, 0,  0,  0, 0, 0, 1,  3,  4, 4, 6, 6,
        6,  6, 5, 3, 1,  0,  0, 0, 0, 0,  2,  3,  5, 5, 6, 6,  6,  5, 4, 2, 1,  0,  0, 0, 0, 0,
        3,  5, 4, 5, 4,  5,  5, 5, 4, 2,  1,  0,  0, 0, 0, 0,  4,  3, 5, 5, 5,  5,  4, 4, 3, 2,
        1,  0, 0, 0, 0,  0,  3, 4, 4, 4,  4,  4,  4, 4, 3, 2,  1,  0, 0, 0, 0,  0};
    std::array<std::uint8_t, width * height> result_mirror{
        0, 0, 0, 0, 0, 0, 0, 0, 4, 8, 12, 12, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 8, 8, 8, 8, 1,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4,  8,  8, 8, 5, 1, 1, 0, 0, 0, 0, 0, 0, 0, 4, 4, 8, 4, 8,
        4, 5, 1, 1, 0, 0, 0, 0, 0, 0, 2,  4,  6, 6, 6, 6, 5, 4, 1, 1, 0, 0, 0, 0, 0, 2, 4, 6, 6,
        4, 4, 4, 5, 4, 1, 1, 0, 0, 0, 0,  2,  4, 4, 4, 4, 6, 6, 3, 3, 3, 1, 1, 0, 0, 0, 2, 4, 4,
        4, 4, 4, 4, 4, 3, 3, 3, 1, 1, 0,  0,  1, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 1, 1, 0, 0, 3,
        4, 4, 4, 4, 3, 2, 3, 4, 4, 4, 3,  1,  1, 0, 0, 4, 3, 2, 3, 4, 4, 4, 4, 3, 2, 3, 3, 1, 1,
        0, 0, 3, 2, 3, 4, 3, 2, 3, 4, 3,  2,  3, 3, 1, 1, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2,
        3, 1, 1, 0, 0, 3, 3, 3, 2, 2, 2,  3,  3, 3, 3, 2, 3, 1, 1, 0, 0, 2, 3, 3, 3, 3, 3, 3, 2,
        2, 2, 2, 3, 1, 1, 0, 0, 3, 3, 2,  2,  2, 3, 3, 2, 2, 2, 2, 3, 1, 1, 0, 0};
    image.print();
    EXPECT_TRUE(arrays_match(result, image.image()));
    EXPECT_TRUE(arrays_match(result_mirror, image.mirror()));
}
