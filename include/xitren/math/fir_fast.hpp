#pragma once

#include <xitren/circular_buffer.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>

namespace xitren::math {

template <std::size_t Order>
class filter : public containers::circular_buffer<std::uint32_t, Order> {
    constexpr static std::uint32_t power_  = 20;
    constexpr static double        factor_ = static_cast<double>(1 << power_);

    using containers::circular_buffer<std::uint32_t, Order>::begin;
    using containers::circular_buffer<std::uint32_t, Order>::end;
    using containers::circular_buffer<std::uint32_t, Order>::full;

public:
    explicit filter(const std::array<double, Order>& table_data) : table_{}
    {
        populate(table_data);
    }

    filter(const std::array<double, Order>&        table_data,
           const std::array<std::uint32_t, Order>& data)
        : table_{}
    {
        populate(table_data);
        (*this) << data;
    }

    filter(const std::array<double, Order>&         table_data,
           const std::array<std::uint32_t, Order>&& data)
        : table_{}
    {
        populate(table_data);
        (*this) << data;
    }

    template <std::size_t Size>
    filter<Size>&
    operator*(const filter<Size>& other)
    {
        static_assert(Order == Size, "Filter lengths do not match!");
        auto other_ptr = other.table_.begin();
        for (auto& item : table_) {
            const double a = (static_cast<double>(item)) / factor_;
            const double b = (static_cast<double>(*other_ptr)) / factor_;
            item           = static_cast<std::uint32_t>((a * b) * factor_);
            other_ptr++;
        }
        return *this;
    }

    template <std::size_t Size>
    filter<Size>&
    operator+(const filter<Size>& other)
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
    operator-(const filter<Size>& other)
    {
        static_assert(Order == Size, "Filter lengths do not match!");
        auto other_ptr = other.table_.begin();
        for (auto& item : table_) {
            item -= *other_ptr;
            other_ptr++;
        }
        return *this;
    }

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

    void
    reset()
    {
        containers::circular_buffer<std::uint32_t, Order>::clear();
    }

    std::array<std::uint32_t, Order>
    table() const
    {
        return table_;
    }

private:
    std::array<std::uint32_t, Order> table_;

    void
    populate(const std::array<double, Order>& table_data)
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
    constexpr static std::array<double, Order + 1>
    prepare_table()
    {
        constexpr double      cutoff = static_cast<double>(Cutoff) / static_cast<double>(Sampling);
        constexpr double      factor = 2.0 * cutoff;
        constexpr std::size_t half   = Order >> 1;
        std::array<double, Order + 1> array{};
        std::size_t                   i = 0;
        for (double& item : array) {
            item = factor * filter<Order + 1>::sinc(factor * ((double)(i++) - (double)half));
        }
        return array;
    }

    lowpass() : filter<Order + 1>{prepare_table()} {}

    template <std::size_t N>
    explicit lowpass(const std::array<std::uint32_t, N>& data)
        : filter<Order + 1>{prepare_table(), data}
    {}

    template <std::size_t N>
    explicit lowpass(const std::array<std::uint32_t, N>&& data)
        : filter<Order + 1>{prepare_table(), data}
    {}

    [[nodiscard]] std::size_t
    order() const
    {
        return Order;
    }
};

template <std::size_t Order, std::size_t Cutoff, std::size_t Sampling>
class highpass : public filter<Order + 1> {
public:
    constexpr static std::array<double, Order + 1>
    prepare_table()
    {
        constexpr double      cutoff = static_cast<double>(Cutoff) / static_cast<double>(Sampling);
        constexpr double      factor = 2.0 * cutoff;
        constexpr std::size_t half   = Order >> 1;
        std::array<double, Order + 1> array{};
        std::size_t                   i = 0;
        for (double& item : array) {
            item = (i == half ? 1.0 : 0.0)
                   - factor * filter<Order + 1>::sinc(factor * ((double)(i++) - (double)half));
        }
        return array;
    }

    highpass() : filter<Order + 1>{prepare_table()} {}

    template <std::size_t N>
    explicit highpass(const std::array<std::uint32_t, N>& data)
        : filter<Order + 1>{prepare_table(), data}
    {}

    template <std::size_t N>
    explicit highpass(const std::array<std::uint32_t, N>&& data)
        : filter<Order + 1>{prepare_table(), data}
    {}

    [[nodiscard]] std::size_t
    order() const
    {
        return Order;
    }
};

template <std::size_t Order, std::size_t LowerCutoff, std::size_t HigherCutoff,
          std::size_t Sampling>
class bandstop : public filter<Order + 1> {
public:
    constexpr static std::array<double, Order + 1>
    prepare_table()
    {
        constexpr auto low  = lowpass<Order, LowerCutoff, Sampling>::prepare_table();
        constexpr auto high = highpass<Order, HigherCutoff, Sampling>::prepare_table();
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

    bandstop() : filter<Order + 1>{prepare_table()} {}

    template <std::size_t N>
    explicit bandstop(const std::array<std::uint32_t, N>& data)
        : filter<Order + 1>{prepare_table(), data}
    {}

    template <std::size_t N>
    explicit bandstop(const std::array<std::uint32_t, N>&& data)
        : filter<Order + 1>{prepare_table(), data}
    {}

    [[nodiscard]] std::size_t
    order() const
    {
        return Order;
    }
};

template <std::size_t Order, std::size_t LowerCutoff, std::size_t HigherCutoff,
          std::size_t Sampling>
class bandpass : public filter<Order + 1> {
public:
    constexpr static std::array<double, Order + 1>
    prepare_table()
    {
        constexpr auto fir = bandstop<Order, LowerCutoff, HigherCutoff, Sampling>::prepare_table();
        constexpr std::size_t         half = Order >> 1;
        std::array<double, Order + 1> array{};
        std::copy(fir.begin(), fir.end(), array.begin());
        std::size_t i = 0;
        for (double& item : array) {
            item = (i++ == half ? 1.0 : 0.0) - item;
        }
        return array;
    }

    bandpass() : filter<Order + 1>{prepare_table()} {}

    template <std::size_t N>
    explicit bandpass(const std::array<std::uint32_t, N>& data)
        : filter<Order + 1>{prepare_table(), data}
    {}

    template <std::size_t N>
    explicit bandpass(const std::array<std::uint32_t, N>&& data)
        : filter<Order + 1>{prepare_table(), data}
    {}

    [[nodiscard]] std::size_t
    order() const
    {
        return Order;
    }
};

template <std::size_t Order>
class moving_average : public filter<Order> {
    constexpr static std::array<double, Order>
    prepare_table()
    {
        std::array<double, Order> array{};
        for (double& item : array) {
            item = (1. / (Order));
        }
        return array;
    }

public:
    moving_average() : filter<Order>{prepare_table()} {}

    template <std::size_t N>
    explicit moving_average(const std::array<std::uint32_t, N>& data)
        : filter<Order>{prepare_table(), data}
    {}

    template <std::size_t N>
    explicit moving_average(const std::array<std::uint32_t, N>&& data)
        : filter<Order>{prepare_table(), data}
    {}

    static std::size_t
    size()
    {
        return Order;
    }
};
}    // namespace loveka::components::math::fir::fast
