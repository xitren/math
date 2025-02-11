#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <vector>

namespace xitren::math {

template <class Type, std::size_t Size>
class matrix : public std::array<std::array<Type, Size>, Size> {

public:
    using data_type = std::array<std::array<Type, Size>, Size>;

    matrix() = default;
    matrix(data_type const& data) : data_type{data} {}

    matrix
    operator*(matrix const& other) const
    {
        data_type ret;
        for (std::size_t i = 0; i < Size; i++) {
            for (std::size_t j = 0; j < Size; j++) {
                for (std::size_t k = 0; k < Size; k++) {
                    ret[i][j] += data_type::operator[](i)[k] * other[k][j];
                }
            }
        }
        return matrix{ret};
    }

    matrix
    operator+(matrix const& other) const
    {
        data_type ret;
        for (std::size_t i = 0; i < Size; i++) {
            for (std::size_t j = 0; j < Size; j++) {
                ret[i][j] = data_type::operator[](i)[j] + other[i][j];
            }
        }
        return matrix{ret};
    }

    matrix
    operator-(matrix const& other) const
    {
        data_type ret;
        for (std::size_t i = 0; i < Size; i++) {
            for (std::size_t j = 0; j < Size; j++) {
                ret[i][j] = data_type::operator[](i)[j] - other[i][j];
            }
        }
        return matrix{ret};
    }

    static matrix
    get_rand_matrix()
    {
        std::srand(std::time({}));    // use current time as seed for random generator
        data_type data_rand;
        for (auto it{data_rand.begin()}; it != data_rand.end(); it++) {
            for (auto it2{(*it).begin()}; it2 != ((*it).end()); it2++) {
                (*it2) = static_cast<Type>(std::rand());
            }
        }
        return matrix{data_rand};
    }
};

}    // namespace xitren::math
