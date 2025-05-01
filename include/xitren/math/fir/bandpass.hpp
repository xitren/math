/*!
_ _
__ _(_) |_ _ _ ___ _ _
\ \ / |  _| '_/ -_) ' \
/_\_\_|\__|_| \___|_||_|
* @date 20.01.2025
*/
#pragma once

#include <xitren/math/fir/bandstop.hpp>
#include <xitren/math/fir/filter.hpp>

#include <array>
#include <cmath>

namespace xitren::math::fir {

template <typename Type, std::size_t Order, std::size_t LowerCutoff, std::size_t HigherCutoff, std::size_t Sampling,
          typename TableType = double>
class bandpass : public filter<Type, Order + 1, TableType> {
    using filter_type = filter<Type, Order + 1, TableType>;

public:
    /**
     * Creates a table of coefficients for a bandpass FIR filter with the given cutoff frequencies and sampling rate.
     * @param cutoff the cutoff frequency of the lowpass section, in samples per second
     * @param sampling_rate the sampling rate of the filter, in samples per second
     * @return the filter table data
     */
    static constexpr filter_type::table_type
    prepare_table()
    {
        constexpr auto        fir  = bandstop<Type, Order, LowerCutoff, HigherCutoff, Sampling>::prepare_table();
        constexpr std::size_t half = Order >> 1;
        typename filter_type::table_type array{};
        std::copy(fir.begin(), fir.end(), array.begin());
        std::size_t i = 0;
        for (double& item : array) {
            item = (i++ == half ? 1.0 : 0.0) - item;
        }
        return array;
    }

    /**
     * Constructs a bandpass FIR filter with the given cutoff frequencies and sampling rate.
     */
    bandpass() : filter_type{prepare_table()} {}

    /**
     * Constructs a bandpass FIR filter with the given cutoff frequencies and sampling rate, and applies it to the given
     * data.
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit bandpass(std::array<double, N> const& data) : filter_type{prepare_table(), data}
    {}

    /**
     * Constructs a bandpass FIR filter with the given cutoff frequencies and sampling rate, and applies it to the given
     * data.
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit bandpass(std::array<double, N> const&& data) : filter_type{prepare_table(), data}
    {}

    /**
     * Returns the order of the bandpass filter.
     * @return the order of the bandpass filter
     */
    [[nodiscard]] std::size_t
    order() const
    {
        return Order;
    }
};
}    // namespace xitren::math::fir
