/*!
_ _
__ _(_) |_ _ _ ___ _ _
\ \ / |  _| '_/ -_) ' \
/_\_\_|\__|_| \___|_||_|
* @date 20.01.2025
*/
#pragma once

#include <xitren/math/fir/filter.hpp>

#include <array>
#include <cmath>

namespace xitren::math::fir {

template <typename Type, std::size_t Order, typename TableType = double>
class moving_average : public filter<Type, Order, TableType> {
    using filter_type = filter<Type, Order, TableType>;

    /**
     * Prepares the filter table data
     * @return the filter table data
     */
    static constexpr filter_type::table_type
    prepare_table()
    {
        typename filter_type::table_type array{};
        for (auto& item : array) {
            item = (1. / static_cast<double>((Order)));
        }
        return array;
    }

public:
    /**
     * Constructs a moving average filter with the given order
     */
    moving_average() : filter_type{prepare_table()} {}

    /**
     * Constructs a moving average filter with the given order and applies it to the given data
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit moving_average(std::array<Type, N> const& data) : filter_type{prepare_table(), data}
    {}

    /**
     * Constructs a moving average filter with the given order and applies it to the given data
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit moving_average(std::array<Type, N> const&& data) : filter_type{prepare_table(), data}
    {}

    /**
     * Returns the order of the moving average filter
     * @return the order of the moving average filter
     */
    static std::size_t
    size()
    {
        return Order;
    }
};
}    // namespace xitren::math::fir
