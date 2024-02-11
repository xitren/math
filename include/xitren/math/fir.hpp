#pragma once

#include <xitren/circular_buffer.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>

namespace xitren::math {

template <std::size_t Order>
class filter : public containers::circular_buffer<double, Order> {
    using containers::circular_buffer<double, Order>::begin;
    using containers::circular_buffer<double, Order>::end;
    using containers::circular_buffer<double, Order>::full;

protected:
    template <typename Real>
    constexpr static Real
    sin_cfrac(Real x2, int k = 2, int n = 40)
    {
        return (n == 0) ? k * (k + 1) - x2
                        : k * (k + 1) - x2 + (k * (k + 1) * x2) / sin_cfrac(x2, k + 2, n - 1);
    }

    template <typename Real>
    constexpr static Real
    wrap(Real x)
    {
        // standardize the angle so that -pi <= x < pi
        return (x <= -M_PI) ? wrap(x + 2 * M_PI) : (x > M_PI) ? wrap(x - 2 * M_PI) : (true) ? x : 0;
    }

    template <typename Real>
    constexpr static Real
    sqr(Real x)
    {
        return x * x;
    }

    template <typename Real>
    constexpr static Real
    sin(Real x)
    {
        return wrap(x) / (1 + sqr(wrap(x)) / sin_cfrac(sqr(wrap(x))));
    }

    template <typename Real>
    constexpr static Real
    sinc(const Real x)
    {
        if (x != 0) {
            const double xpi = M_PI * x;
            return sin(xpi) / xpi;
        }
        return 1.0;
    }

    constexpr static double
    i0(const double x)
    {
        double       f  = 1;
        const double x2 = x * x * 0.25;
        double       xc = x2;
        double       v  = 1 + x2;
        for (int i = 2; i < 100; i++) {
            f *= i;
            xc *= x2;
            const double a = xc / (f * f);
            v += a;
            if (a < 1e-20)
                break;
        }
        return v;
    }

    template <std::size_t Size>
    constexpr static void
    normalize(std::array<double, Size>& win)
    {
        double sum = 0;
        for (auto& item : win) {
            sum += item;
        }
        if (sum != 0) {
            for (auto& item : win) {
                item /= sum;
            }
        }
    }

    template <std::size_t Size>
    constexpr static void
    window_kaiser(std::array<double, Size>& win, const double transitionWidth,
                  const double attenuation, const double fs)
    {
        const double tw = 2.0 * M_PI * transitionWidth / fs;
        std::int32_t m  = (attenuation <= 21) ? ((int)::ceil(5.79 / tw))
                                              : ((int)::ceil((attenuation - 7.95) / (2.285 * tw)));
        if ((m & 1) == 0)
            m++;
        double beta
            = (attenuation <= 21)
                  ? (0)
                  : ((attenuation <= 50)
                         ? (0.5842 * ::pow(attenuation - 21, 0.4) + 0.07886 * (attenuation - 21))
                         : (0.1102 * (attenuation - 8.7)));
        const double i0b = i0(beta);
        for (int n = 0; n < m; n++) {
            const double v = beta * ::sqrt(1.0 - ::pow(2.0 * n / (m - 1) - 1.0, 2));
            win[n]         = i0(v) / i0b;
        }
    }

    template <std::size_t Size>
    constexpr static void
    window_blackman(std::array<double, Size>& win)
    {
        for (int i = 0; i < Size; i++) {
            win[i] *= 0.42 - 0.5 * ::cos(2.0 * M_PI * i / (Size - 1))
                      + 0.08 * ::cos(4.0 * M_PI * i / (Size - 1));
        }
    }

    template <std::size_t Size>
    constexpr static void
    window_sinc(std::array<double, Size>& win)
    {
        const std::size_t m = Size - 1;
        std::size_t       i = 0;
        for (auto& item : win) {
            item *= sinc(2.0 * (double)(i++) / m - 1.0);
        }
    }

    template <std::size_t Size>
    constexpr static void
    window_hanning(std::array<double, Size>& win)
    {
        const std::size_t m = Size - 1;
        std::size_t       i = 0;
        for (auto& item : win) {
            item *= 0.5 - 0.5 * ::cos(2.0 * M_PI * (double)(i++) / m);
        }
    }

    template <std::size_t Size>
    constexpr static void
    window_hamming(std::array<double, Size>& win)
    {
        const std::size_t m = Size - 1;
        std::size_t       i = 0;
        for (auto& item : win) {
            item *= 0.54 - 0.46 * ::cos(2.0 * M_PI * (double)(i++) / m);
        }
    }

public:
    constexpr explicit filter(const std::array<double, Order>& table_data) : table_{table_data} {}

    constexpr filter(const std::array<double, Order>& table_data,
                     const std::array<double, Order>& data)
        : table_{table_data}
    {
        (*this) << data;
    }

    filter(const std::array<double, Order>& table_data, const std::array<double, Order>&& data)
        : table_{table_data}
    {
        (*this) << data;
    }

    template <std::size_t Size>
    filter<Size>&
    operator*(const filter<Size>& other)
    {
        static_assert(Order == Size, "Filter lengths do not match!");
        auto other_ptr = other.table_.begin();
        for (auto& item : table_) {
            item *= *other_ptr;
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

    double
    value(double val)
    {
        containers::circular_buffer<double, Order>::push(val);
        if (!full())
            return 0.;
        auto   it      = begin();
        double ret_val = 0.;
        for (auto& item : table_) {
            ret_val += item * (*it);
            it++;
        }
        return ret_val;
    }

    void
    reset()
    {
        containers::circular_buffer<double, Order>::clear();
    }

    std::array<double, Order>
    table() const
    {
        return table_;
    }

private:
    std::array<double, Order> table_;
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

    constexpr lowpass() : filter<Order + 1>{prepare_table()} {}

    template <std::size_t N>
    constexpr explicit lowpass(const std::array<double, N>& data)
        : filter<Order + 1>{prepare_table(), data}
    {}

    template <std::size_t N>
    constexpr explicit lowpass(const std::array<double, N>&& data)
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
                   - factor * filter<Order + 1>::sinc(factor * ((double)(i) - (double)half));
            i++;
        }
        return array;
    }

    highpass() : filter<Order + 1>{prepare_table()} {}

    template <std::size_t N>
    explicit highpass(const std::array<double, N>& data) : filter<Order + 1>{prepare_table(), data}
    {}

    template <std::size_t N>
    explicit highpass(const std::array<double, N>&& data) : filter<Order + 1>{prepare_table(), data}
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
        constexpr std::array<double, Order + 1> low
            = lowpass<Order, LowerCutoff, Sampling>::prepare_table();
        constexpr std::array<double, Order + 1> high
            = highpass<Order, HigherCutoff, Sampling>::prepare_table();
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
    explicit bandstop(const std::array<double, N>& data) : filter<Order + 1>{prepare_table(), data}
    {}

    template <std::size_t N>
    explicit bandstop(const std::array<double, N>&& data) : filter<Order + 1>{prepare_table(), data}
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
    explicit bandpass(const std::array<double, N>& data) : filter<Order + 1>{prepare_table(), data}
    {}

    template <std::size_t N>
    explicit bandpass(const std::array<double, N>&& data) : filter<Order + 1>{prepare_table(), data}
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
    explicit moving_average(const std::array<double, N>& data)
        : filter<Order>{prepare_table(), data}
    {}

    template <std::size_t N>
    explicit moving_average(const std::array<double, N>&& data)
        : filter<Order>{prepare_table(), data}
    {}

    static std::size_t
    size()
    {
        return Order;
    }
};
}    // namespace loveka::components::math::fir
