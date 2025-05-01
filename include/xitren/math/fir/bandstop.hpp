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

template <typename Type, std::size_t Order, std::size_t LowerCutoff, std::size_t HigherCutoff, std::size_t Sampling,
          typename TableType = double>
class bandstop : public filter<Type, Order + 1, TableType> {
    using filter_type = filter<Type, Order + 1, TableType>;

public:
    /**
     * Creates a table of coefficients for a bandstop FIR filter with the given cutoff frequencies and sampling rate.
     * @param cutoff the cutoff frequency of the lowpass section, in samples per second
     * @param sampling_rate the sampling rate of the filter, in samples per second
     * @return the filter table data
     */
    static constexpr filter_type::table_type
    prepare_table()
    {
        constexpr typename filter_type::table_type low
            = lowpass<Type, Order, LowerCutoff, Sampling, TableType>::prepare_table();
        constexpr typename filter_type::table_type high
            = highpass<Type, Order, HigherCutoff, Sampling, TableType>::prepare_table();
        std::array<double, Order + 1> array{};
        auto                          low_ptr = low.begin();
        for (double& item : array) {
            item = *low_ptr;
            ++low_ptr;
        }
        auto high_ptr = high.begin();
        for (double& item : array) {
            item += *high_ptr;
            ++high_ptr;
        }
        return array;
    }

    /**
     * Constructs a bandstop FIR filter with the given cutoff frequencies and sampling rate.
     */
    bandstop() : filter_type{prepare_table()} {}

    /**
     * Constructs a bandstop FIR filter with the given cutoff frequencies and sampling rate, and applies it to the given
     * data.
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit bandstop(std::array<double, N> const& data) : filter_type{prepare_table(), data}
    {}

    /**
     * Constructs a bandstop FIR filter with the given cutoff frequencies and sampling rate, and applies it to the given
     * data.
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit bandstop(std::array<double, N> const&& data) : filter_type{prepare_table(), data}
    {}

    /**
     * Returns the order of the bandstop filter.
     * @return the order of the bandstop filter
     */
    [[nodiscard]] std::size_t
    order() const
    {
        return Order;
    }
};
}    // namespace xitren::math::fir
