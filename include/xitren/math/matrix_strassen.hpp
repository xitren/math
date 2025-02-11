#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <vector>

namespace xitren::math {

template <class Type, std::size_t Size>
class matrix_strassen {
    static_assert((Size & (Size - 1)) == 0, "Should be power of 2!");

public:
    using quarter_type      = matrix_strassen<Type, Size / 2>;
    using data_type         = std::array<Type, Size * Size>;
    using quarter_data_type = std::array<Type, Size * Size / 4>;

    matrix_strassen() = default;
    matrix_strassen(data_type const& data)
        : a_{get_init(data, 0)}, b_{get_init(data, 1)}, c_{get_init(data, 2)}, d_{get_init(data, 3)}
    {}

    auto
    get(std::size_t row, std::size_t column) const
    {
        if (row < (Size / 2)) {
            if (column < (Size / 2)) {
                return a_.get(row % (Size / 2), column % (Size / 2));
            } else {
                return b_.get(row % (Size / 2), column % (Size / 2));
            }
        } else {
            if (column < (Size / 2)) {
                return c_.get(row % (Size / 2), column % (Size / 2));
            } else {
                return d_.get(row % (Size / 2), column % (Size / 2));
            }
        }
    }

    matrix_strassen
    operator*(matrix_strassen const& other) const
    {
        auto const H1 = (a_ + d_) * (other.a_ + other.d_);
        auto const H2 = (c_ + d_) * other.a_;
        auto const H3 = a_ * (other.b_ - other.d_);
        auto const H4 = d_ * (other.c_ - other.a_);
        auto const H5 = (a_ + b_) * other.d_;
        auto const H6 = (c_ - a_) * (other.a_ + other.b_);
        auto const H7 = (b_ - d_) * (other.c_ + other.d_);
        return matrix_strassen{H1 + H4 - H5 + H7, H3 + H5, H2 + H4, H1 + H3 - H2 + H6};
    }

    matrix_strassen
    operator+(matrix_strassen const& other) const
    {
        return matrix_strassen{a_ + other.a_, b_ + other.b_, c_ + other.c_, d_ + other.d_};
    }

    matrix_strassen
    operator-(matrix_strassen const& other) const
    {
        return matrix_strassen{a_ + other.a_, b_ + other.b_, c_ + other.c_, d_ + other.d_};
    }

    static matrix_strassen
    get_rand_matrix()
    {
        std::srand(std::time({}));    // use current time as seed for random generator
        data_type data_rand;
        for (auto it{data_rand.begin()}; it != data_rand.end(); it++) {
            (*it) = static_cast<Type>(std::rand());
        }
        return matrix_strassen{data_rand};
    }

private:
    const quarter_type a_{};
    const quarter_type b_{};
    const quarter_type c_{};
    const quarter_type d_{};

    matrix_strassen(quarter_type const& m_a, quarter_type const& m_b, quarter_type const& m_c, quarter_type const& m_d)
        : a_{m_a}, b_{m_b}, c_{m_c}, d_{m_d}
    {}

    static constexpr quarter_data_type
    get_init(data_type const& data, std::size_t k)
    {
        quarter_data_type ret;
        std::size_t       x{};
        std::size_t       y{};
        std::size_t       j{};
        switch (k) {
        case 1:
            x = Size / 2;
            break;
        case 2:
            y = Size / 2;
            break;
        case 3:
            x = Size / 2;
            y = Size / 2;
            break;
        default:
            break;
        }
        for (std::size_t l{}; l < (Size / 2); l++) {
            for (std::size_t m{}; m < (Size / 2); m++) {
                ret[j++] = data[(l + y) * Size + m + x];
            }
        }
        return ret;
    }
};

template <class Type>
class matrix_strassen<Type, 2> : public std::array<Type, 4> {

public:
    using data_type = std::array<Type, 4>;

    matrix_strassen() = default;
    matrix_strassen(const data_type data) : data_type{data} {}

    auto
    get(std::size_t row, std::size_t column) const
    {
        return data_type::operator[](row * 2 + column);
    }

    matrix_strassen
    operator*(matrix_strassen const& other) const
    {
        Type const& a = data_type::operator[](0);
        Type const& b = data_type::operator[](1);
        Type const& c = data_type::operator[](2);
        Type const& d = data_type::operator[](3);
        Type const& A = other[0];
        Type const& C = other[1];
        Type const& B = other[2];
        Type const& D = other[3];

        const Type t = a * A;
        const Type u = (c - a) * (C - D);
        const Type v = (c + d) * (C - A);
        const Type w = t + (c + d - a) * (A + D - C);

        return matrix_strassen{typename matrix_strassen::data_type{
            {t + b * B, w + v + (a + b - c - d) * D, w + u + d * (B + C - A - D), w + u + v}}};
    }

    matrix_strassen
    operator+(matrix_strassen const& other) const
    {
        return matrix_strassen{typename matrix_strassen::data_type{
            {data_type::operator[](0) + other[0], data_type::operator[](1) + other[1],
             data_type::operator[](2) + other[2], data_type::operator[](3) + other[3]}}};
    }

    matrix_strassen
    operator-(matrix_strassen const& other) const
    {
        return matrix_strassen{typename matrix_strassen::data_type{
            {data_type::operator[](0) - other[0], data_type::operator[](1) - other[1],
             data_type::operator[](2) - other[2], data_type::operator[](3) - other[3]}}};
    }

private:
};

}    // namespace xitren::math
