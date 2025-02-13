#pragma once

#include <xitren/math/matrix_strassen.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <memory>
#include <vector>

namespace xitren::math {

template <class Type, std::size_t Rows, std::size_t Columns>
class matrix {

    static constexpr std::size_t
    get_max_item()
    {
        std::array<std::size_t, 7> sizes{{128, 64, 32, 16, 8, 4, 2}};

        std::size_t max_row{};
        for (auto& size : sizes) {
            if ((Rows / size) != 0) {
                max_row = size;
            }
            break;
        }
        std::size_t max_col{};
        for (auto& size : sizes) {
            if ((Columns / size) != 0) {
                max_col = size;
            }
            break;
        }
        return std::min(max_row, max_col);
    }

public:
    static constexpr auto batch_value = get_max_item();

    using init_type        = std::vector<Type>;
    using batch_type       = matrix_strassen<Type, batch_value>;
    using batch_store_type = std::shared_ptr<batch_type>;

    static constexpr auto batch_rows    = Rows / batch_value;
    static constexpr auto batch_columns = Columns / batch_value;
    static constexpr auto rest_rows     = Rows % batch_value;
    static constexpr auto rest_columns  = Columns % batch_value;

    using data_type         = std::array<std::array<batch_store_type, batch_columns>, Rows / batch_value>;
    using rest_rows_type    = std::array<Type>;
    using rest_columns_type = std::array<Type>;

    matrix() = default;
    matrix(init_type const& data) {}

    static constexpr std::array<Type, batch_type>
    get_batch_init(data_type const& data, std::size_t k)
    {}

    // template <std::size_t ColumnsOther>
    // matrix_classic
    // operator*(matrix_classic<Type, Columns, ColumnsOther> const& other) const
    // {
    //     data_type ret;
    //     for (std::size_t i = 0; i < Rows; i++) {
    //         for (std::size_t j = 0; j < ColumnsOther; j++) {
    //             for (std::size_t k = 0; k < Columns; k++) {
    //                 ret[i][j] += data_type::operator[](i)[k] * other[k][j];
    //             }
    //         }
    //     }
    //     return matrix_classic{ret};
    // }

    // matrix_classic
    // operator+(matrix_classic const& other) const
    // {
    //     data_type ret;
    //     for (std::size_t i = 0; i < Rows; i++) {
    //         for (std::size_t j = 0; j < Columns; j++) {
    //             ret[i][j] = data_type::operator[](i)[j] + other[i][j];
    //         }
    //     }
    //     return matrix_classic{ret};
    // }

    // matrix_classic
    // operator-(matrix_classic const& other) const
    // {
    //     data_type ret;
    //     for (std::size_t i = 0; i < Rows; i++) {
    //         for (std::size_t j = 0; j < Columns; j++) {
    //             ret[i][j] = data_type::operator[](i)[j] - other[i][j];
    //         }
    //     }
    //     return matrix_classic{ret};
    // }

    // static matrix_classic
    // get_rand_matrix()
    // {
    //     std::srand(std::time({}));    // use current time as seed for random generator
    //     data_type data_rand;
    //     for (auto it{data_rand.begin()}; it != data_rand.end(); it++) {
    //         for (auto it2{(*it).begin()}; it2 != ((*it).end()); it2++) {
    //             (*it2) = static_cast<Type>(std::rand());
    //         }
    //     }
    //     return matrix_classic{data_rand};
    // }

private:
    data_type         batched_section;
    rest_rows_type    rest_rows_section;
    rest_columns_type rest_columns_section;
};

}    // namespace xitren::math
