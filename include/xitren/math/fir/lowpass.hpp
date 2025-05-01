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

template <typename Type, std::size_t Order, std::size_t Cutoff, std::size_t Sampling, typename TableType = double>
class lowpass : public filter<Type, Order + 1, TableType> {
    using filter_type = filter<Type, Order + 1, TableType>;

public:
    /**
     * Creates a lowpass filter with the given cutoff frequency and sampling rate
     * @param cutoff the cutoff frequency of the filter, in samples per second
     * @param sampling_rate the sampling rate of the filter, in samples per second
     */
    static constexpr filter_type::table_type
    prepare_table()
    {
        // calculate the cutoff frequency in terms of the filter's sampling rate
        constexpr double cutoff = static_cast<double>(Cutoff) / static_cast<double>(Sampling);
        // calculate the filter's factor, which determines the amount of attenuation
        constexpr double                 factor = 2.0 * cutoff;
        constexpr std::size_t            half   = Order >> 1;
        typename filter_type::table_type array{};
        // loop through each element in the filter's table
        std::size_t i{};
        for (auto& item : array) {
            // calculate the current element of the filter's table using the sinc function
            item = factor * filter_type::sinc(factor * (static_cast<double>(i++) - static_cast<double>(half)));
        }
        return array;
    }

    /**
     * Constructs a lowpass filter with the given cutoff frequency and sampling rate
     */
    constexpr lowpass() : filter_type{prepare_table()} {}

    /**
     * Constructs a lowpass filter with the given cutoff frequency and sampling rate, and applies it to the given data
     * @param data the data to filter
     */
    template <std::size_t N>
    constexpr explicit lowpass(std::array<double, N> const& data) : filter_type{prepare_table(), data}
    {}

    /**
     * Constructs a lowpass filter with the given cutoff frequency and sampling rate, and applies it to the given data
     * @param data the data to filter
     */
    template <std::size_t N>
    constexpr explicit lowpass(std::array<double, N> const&& data) : filter_type{prepare_table(), data}
    {}

    /**
     * Returns the order of the lowpass filter
     * @return the order of the lowpass filter
     */
    [[nodiscard]] std::size_t
    order() const
    {
        return Order;
    }
};

}    // namespace xitren::math::fir
