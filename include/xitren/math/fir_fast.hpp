#pragma once

#include <xitren/circular_buffer.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>

namespace xitren::math {

template <std::size_t Order>
class filter : public containers::circular_buffer<std::uint32_t, Order> {
    static constexpr std::uint32_t power_  = 20;
    static constexpr double        factor_ = static_cast<double>(1 << power_);

    using containers::circular_buffer<std::uint32_t, Order>::begin;
    using containers::circular_buffer<std::uint32_t, Order>::end;
    using containers::circular_buffer<std::uint32_t, Order>::full;

public:
    /**
     * Constructs a filter with the given table data
     * @param table_data the table data to use for the filter
     */
    explicit filter(std::array<double, Order> const& table_data) : table_{} { populate(table_data); }

    /**
     * Constructs a filter with the given table data and data
     * @param table_data the table data to use for the filter
     * @param data the data to filter
     */
    filter(std::array<double, Order> const& table_data, std::array<std::uint32_t, Order> const& data) : table_{}
    {
        populate(table_data);
        (*this) << data;
    }

    /**
     * Constructs a filter with the given table data and rvalue data
     * @param table_data the table data to use for the filter
     * @param data the rvalue data to filter
     */
    filter(std::array<double, Order> const& table_data, std::array<std::uint32_t, Order> const&& data) : table_{}
    {
        populate(table_data);
        (*this) << data;
    }

    template <std::size_t Size>
    filter<Size>&
    operator*(filter<Size> const& other)
    {
        static_assert(Order == Size, "Filter lengths do not match!");
        auto other_ptr = other.table_.begin();
        for (auto& item : table_) {
            double const a = (static_cast<double>(item)) / factor_;
            double const b = (static_cast<double>(*other_ptr)) / factor_;
            item           = static_cast<std::uint32_t>((a * b) * factor_);
            other_ptr++;
        }
        return *this;
    }

    template <std::size_t Size>
    filter<Size>&
    operator+(filter<Size> const& other)
    {
        static_assert(Order == Size, "Filter lengths do not match!");
        auto other_ptr = other.table_.begin();
        for (auto& item : table_) {
            item += *other_ptr;
            other_ptr++;
        }
        return *this;
    }

    template <std::size_t Size>
    filter<Size>&
    operator-(filter<Size> const& other)
    {
        static_assert(Order == Size, "Filter lengths do not match!");
        auto other_ptr = other.table_.begin();
        for (auto& item : table_) {
            item -= *other_ptr;
            other_ptr++;
        }
        return *this;
    }

    /**
     * Applies the filter to a new data point
     * @param val the new data point
     * @return the filtered data point
     */
    std::uint32_t
    value(std::uint32_t val)
    {
        containers::circular_buffer<std::uint32_t, Order>::push(val);
        if (!full())
            return 0;
        auto          it      = begin();
        std::uint64_t ret_val = 0;
        for (auto& item : table_) {
            ret_val += item * (*it);
            it++;
        }
        ret_val = ret_val >> power_;
        return static_cast<std::uint32_t>(ret_val);
    }

    /**
     * Resets the filter state
     */
    void
    reset()
    {
        containers::circular_buffer<std::uint32_t, Order>::clear();
    }

    /**
     * Returns the filter table data
     * @return the filter table data
     */
    std::array<std::uint32_t, Order>
    table() const
    {
        return table_;
    }

private:
    std::array<std::uint32_t, Order> table_;

    void
    populate(std::array<double, Order> const& table_data)
    {
        auto other_ptr = table_data.begin();
        for (auto& item : table_) {
            item = static_cast<std::uint32_t>((*other_ptr) * factor_);
            other_ptr++;
        }
    }
};

template <std::size_t Order, std::size_t Cutoff, std::size_t Sampling>
class lowpass : public filter<Order + 1> {
public:
    /**
     * Creates a lowpass filter with the given cutoff frequency and sampling rate
     * @param cutoff the cutoff frequency of the filter, in samples per second
     * @param sampling_rate the sampling rate of the filter, in samples per second
     */
    static constexpr std::array<double, Order + 1>
    prepare_table()
    {
        // calculate the cutoff frequency in terms of the filter's sampling rate
        constexpr double cutoff = static_cast<double>(Cutoff) / static_cast<double>(Sampling);
        // calculate the filter's factor, which determines the amount of attenuation
        constexpr double              factor = 2.0 * cutoff;
        constexpr std::size_t         half   = Order >> 1;
        std::array<double, Order + 1> array{};
        // loop through each element in the filter's table
        std::size_t i = 0;
        for (double& item : array) {
            // calculate the current element of the filter's table using the sinc function
            item = factor * filter<Order + 1>::sinc(factor * ((double)(i++) - (double)half));
        }
        return array;
    }

    /**
     * Constructs a lowpass filter with the given cutoff frequency and sampling rate
     */
    lowpass() : filter<Order + 1> { prepare_table() }
    {}

    /**
     * Constructs a lowpass filter with the given cutoff frequency and sampling rate, and applies it to the given data
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit lowpass(std::array<std::uint32_t, N> const& data) : filter<Order + 1>
    {
        prepare_table(), data
    }
    {}

    /**
     * Constructs a lowpass filter with the given cutoff frequency and sampling rate, and applies it to the given data
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit lowpass(std::array<std::uint32_t, N> const&& data) : filter<Order + 1>
    {
        prepare_table(), data
    }
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

template <std::size_t Order, std::size_t Cutoff, std::size_t Sampling>
class highpass : public filter<Order + 1> {
public:
    /**
     * Creates a highpass filter with the given cutoff frequency and sampling rate
     * @param cutoff the cutoff frequency of the filter, in samples per second
     * @param sampling_rate the sampling rate of the filter, in samples per second
     */
    static constexpr std::array<double, Order + 1>
    prepare_table()
    {
        // calculate the cutoff frequency in terms of the filter's sampling rate
        constexpr double cutoff = static_cast<double>(Cutoff) / static_cast<double>(Sampling);
        // calculate the filter's factor, which determines the amount of attenuation
        constexpr double              factor = 2.0 * cutoff;
        constexpr std::size_t         half   = Order >> 1;
        std::array<double, Order + 1> array{};
        // loop through each element in the filter's table
        std::size_t i = 0;
        for (double& item : array) {
            // calculate the current element of the filter's table using the sinc function
            item = (i == half ? 1.0 : 0.0) - factor * filter<Order + 1>::sinc(factor * ((double)(i++) - (double)half));
        }
        return array;
    }

    /**
     * Constructs a highpass filter with the given cutoff frequency and sampling rate
     */
    highpass() : filter<Order + 1> { prepare_table() }
    {}

    /**
     * Constructs a highpass filter with the given cutoff frequency and sampling rate, and applies it to the given data
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit highpass(std::array<std::uint32_t, N> const& data) : filter<Order + 1>
    {
        prepare_table(), data
    }
    {}

    /**
     * Constructs a highpass filter with the given cutoff frequency and sampling rate, and applies it to the given data
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit highpass(std::array<std::uint32_t, N> const&& data) : filter<Order + 1>
    {
        prepare_table(), data
    }
    {}

    /**
     * Returns the order of the highpass filter
     * @return the order of the highpass filter
     */
    [[nodiscard]] std::size_t
    order() const
    {
        return Order;
    }
};

template <std::size_t Order, std::size_t LowerCutoff, std::size_t HigherCutoff, std::size_t Sampling>
class bandstop : public filter<Order + 1> {
public:
    /**
     * Creates a table of coefficients for a bandstop FIR filter with the given cutoff frequencies and sampling rate.
     * @param cutoff the cutoff frequency of the lowpass section, in samples per second
     * @param sampling_rate the sampling rate of the filter, in samples per second
     * @return the filter table data
     */
    static constexpr std::array<double, Order + 1>
    prepare_table()
    {
        constexpr auto                low  = lowpass<Order, LowerCutoff, Sampling>::prepare_table();
        constexpr auto                high = highpass<Order, HigherCutoff, Sampling>::prepare_table();
        std::array<double, Order + 1> array{};
        auto                          low_ptr = low.begin();
        for (double& item : array) {
            item = *low_ptr;
            low_ptr++;
        }
        auto high_ptr = high.begin();
        for (double& item : array) {
            item += *high_ptr;
            high_ptr++;
        }
        return array;
    }

    /**
     * Constructs a bandstop FIR filter with the given cutoff frequencies and sampling rate.
     */
    bandstop() : filter<Order + 1> { prepare_table() }
    {}

    /**
     * Constructs a bandstop FIR filter with the given cutoff frequencies and sampling rate, and applies it to the given
     * data.
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit bandstop(std::array<std::uint32_t, N> const& data) : filter<Order + 1>
    {
        prepare_table(), data
    }
    {}

    /**
     * Constructs a bandstop FIR filter with the given cutoff frequencies and sampling rate, and applies it to the given
     * data.
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit bandstop(std::array<std::uint32_t, N> const&& data) : filter<Order + 1>
    {
        prepare_table(), data
    }
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

template <std::size_t Order, std::size_t LowerCutoff, std::size_t HigherCutoff, std::size_t Sampling>
class bandpass : public filter<Order + 1> {
public:
    /**
     * Creates a table of coefficients for a bandpass FIR filter with the given cutoff frequencies and sampling rate.
     * @param cutoff the cutoff frequency of the lowpass section, in samples per second
     * @param sampling_rate the sampling rate of the filter, in samples per second
     * @return the filter table data
     */
    static constexpr std::array<double, Order + 1>
    prepare_table()
    {
        constexpr auto                fir  = bandstop<Order, LowerCutoff, HigherCutoff, Sampling>::prepare_table();
        constexpr std::size_t         half = Order >> 1;
        std::array<double, Order + 1> array{};
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
    bandpass() : filter<Order + 1> { prepare_table() }
    {}

    /**
     * Constructs a bandpass FIR filter with the given cutoff frequencies and sampling rate, and applies it to the given
     * data.
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit bandpass(std::array<std::uint32_t, N> const& data) : filter<Order + 1>
    {
        prepare_table(), data
    }
    {}

    /**
     * Constructs a bandpass FIR filter with the given cutoff frequencies and sampling rate, and applies it to the given
     * data.
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit bandpass(std::array<std::uint32_t, N> const&& data) : filter<Order + 1>
    {
        prepare_table(), data
    }
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

template <std::size_t Order>
class moving_average : public filter<Order> {
    /**
     * Prepares the filter table data
     * @return the filter table data
     */
    static constexpr std::array<double, Order>
    prepare_table()
    {
        std::array<double, Order> array{};
        for (double& item : array) {
            item = (1. / (Order));
        }
        return array;
    }

public:
    /**
     * Constructs a moving average filter with the given order
     */
    moving_average() : filter<Order> { prepare_table() }
    {}

    /**
     * Constructs a moving average filter with the given order and applies it to the given data
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit moving_average(std::array<std::uint32_t, N> const& data) : filter<Order>
    {
        prepare_table(), data
    }
    {}

    /**
     * Constructs a moving average filter with the given order and applies it to the given data
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit moving_average(std::array<std::uint32_t, N> const&& data) : filter<Order>
    {
        prepare_table(), data
    }
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
}    // namespace xitren::math
