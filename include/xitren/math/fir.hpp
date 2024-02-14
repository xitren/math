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
    /**
     * Calculates the k-th order modified Bessel function of the first kind
     * @param x2 the argument of the function
     * @param k the order of the function
     * @param n the number of terms used in the Taylor series approximation
     * @return the k-th order modified Bessel function of the first kind evaluated at x
     */
    template <typename Real>
    static constexpr Real
    sin_cfrac(Real x2, int k = 2, int n = 40)
    {
        return (n == 0) ? k * (k + 1) - x2 : k * (k + 1) - x2 + (k * (k + 1) * x2) / sin_cfrac(x2, k + 2, n - 1);
    }

    /**
     * Wraps the given angle x so that it lies in the range [-pi, pi)
     * @param x the angle to wrap
     * @return the wrapped angle
     */
    template <typename Real>
    static constexpr Real
    wrap(Real x)
    {
        // standardize the angle so that -pi <= x < pi
        return (x <= -M_PI) ? wrap(x + 2 * M_PI) : (x > M_PI) ? wrap(x - 2 * M_PI) : (true) ? x : 0;
    }

    /**
     * Calculates the square of the given value
     * @param x the value to square
     * @return the square of x
     */
    template <typename Real>
    static constexpr Real
    sqr(Real x)
    {
        return x * x;
    }

    /**
     * Calculates the sine of the given angle
     * @param x the angle to calculate the sine of
     * @return the sine of x
     */
    template <typename Real>
    static constexpr Real
    sin(Real x)
    {
        return wrap(x) / (1 + sqr(wrap(x)) / sin_cfrac(sqr(wrap(x))));
    }

    /**
     * Calculates the sinc function of the given argument
     * @param x the argument of the sinc function
     * @return the sinc function of x
     */
    template <typename Real>
    static constexpr Real
    sinc(Real const x)
    {
        if (x != 0) {
            double const xpi = M_PI * x;
            return sin(xpi) / xpi;
        }
        return 1.0;
    }

    /**
     * Calculates the modified Bessel function of the first kind of order 0
     * @param x the argument of the function
     * @return the modified Bessel function of the first kind of order 0 evaluated at x
     */
    static constexpr double
    i0(double const x)
    {
        double       f  = 1;
        double const x2 = x * x * 0.25;
        double       xc = x2;
        double       v  = 1 + x2;
        for (int i = 2; i < 100; i++) {
            f *= i;
            xc *= x2;
            double const a = xc / (f * f);
            v += a;
            if (a < 1e-20)
                break;
        }
        return v;
    }

    /**
     * Normalizes the given window so that its sum is 1
     * @param win the window to normalize
     */
    template <std::size_t Size>
    static constexpr void
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

    /**
     * Creates a Kaiser window with the given parameters
     * @param win the window to create
     * @param transitionWidth the transition width of the window
     * @param attenuation the attenuation of the window
     * @param fs the sampling frequency of the window
     */
    template <std::size_t Size>
    static constexpr void
    window_kaiser(std::array<double, Size>& win, double const transitionWidth, double const attenuation,
                  double const fs)
    {
        double const tw = 2.0 * M_PI * transitionWidth / fs;
        std::int32_t m
            = (attenuation <= 21) ? ((int)::ceil(5.79 / tw)) : ((int)::ceil((attenuation - 7.95) / (2.285 * tw)));
        if ((m & 1) == 0)
            m++;
        double beta
            = (attenuation <= 21)
                  ? (0)
                  : ((attenuation <= 50) ? (0.5842 * ::pow(attenuation - 21, 0.4) + 0.07886 * (attenuation - 21))
                                         : (0.1102 * (attenuation - 8.7)));
        double const i0b = i0(beta);
        for (int n = 0; n < m; n++) {
            double const v = beta * ::sqrt(1.0 - ::pow(2.0 * n / (m - 1) - 1.0, 2));
            win[n]         = i0(v) / i0b;
        }
    }

    /**
     * Creates a Blackman window with the given parameters
     * @param win the window to create
     */
    template <std::size_t Size>
    static constexpr void
    window_blackman(std::array<double, Size>& win)
    {
        for (int i = 0; i < Size; i++) {
            win[i] *= 0.42 - 0.5 * ::cos(2.0 * M_PI * i / (Size - 1)) + 0.08 * ::cos(4.0 * M_PI * i / (Size - 1));
        }
    }

    /**
     * Creates a Sinc window with the given parameters
     * @param win the window to create
     */
    template <std::size_t Size>
    static constexpr void
    window_sinc(std::array<double, Size>& win)
    {
        std::size_t const m = Size - 1;
        std::size_t       i = 0;
        for (auto& item : win) {
            item *= sinc(2.0 * (double)(i++) / m - 1.0);
        }
    }

    /**
     * Creates a Hanning window with the given parameters
     * @param win the window to create
     */
    template <std::size_t Size>
    static constexpr void
    window_hanning(std::array<double, Size>& win)
    {
        std::size_t const m = Size - 1;
        std::size_t       i = 0;
        for (auto& item : win) {
            item *= 0.5 - 0.5 * ::cos(2.0 * M_PI * (double)(i++) / m);
        }
    }

    /**
     * Creates a Hamming window with the given parameters
     * @param win the window to create
     */
    template <std::size_t Size>
    static constexpr void
    window_hamming(std::array<double, Size>& win)
    {
        std::size_t const m = Size - 1;
        std::size_t       i = 0;
        for (auto& item : win) {
            item *= 0.54 - 0.46 * ::cos(2.0 * M_PI * (double)(i++) / m);
        }
    }

public:
    /**
     * Constructs a filter with the given table data
     * @param table_data the table data to use for the filter
     */
    constexpr explicit filter(std::array<double, Order> const& table_data) : table_{table_data} {}

    /**
     * Constructs a filter with the given table data and data
     * @param table_data the table data to use for the filter
     * @param data the data to filter
     */
    constexpr filter(std::array<double, Order> const& table_data, std::array<double, Order> const& data)
        : table_{table_data}
    {
        (*this) << data;
    }

    /**
     * Constructs a filter with the given table data and rvalue data
     * @param table_data the table data to use for the filter
     * @param data the rvalue data to filter
     */
    filter(std::array<double, Order> const& table_data, std::array<double, Order> const&& data) : table_{table_data}
    {
        (*this) << data;
    }

    template <std::size_t Size>
    filter<Size>&
    operator*(filter<Size> const& other)
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

    /**
     * Resets the filter state
     */
    void
    reset()
    {
        containers::circular_buffer<double, Order>::clear();
    }

    /**
     * Returns the filter table data
     * @return the filter table data
     */
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
    constexpr lowpass() : filter<Order + 1> { prepare_table() }
    {}

    /**
     * Constructs a lowpass filter with the given cutoff frequency and sampling rate, and applies it to the given data
     * @param data the data to filter
     */
    template <std::size_t N>
    constexpr explicit lowpass(std::array<double, N> const& data) : filter<Order + 1>
    {
        prepare_table(), data
    }
    {}

    /**
     * Constructs a lowpass filter with the given cutoff frequency and sampling rate, and applies it to the given data
     * @param data the data to filter
     */
    template <std::size_t N>
    constexpr explicit lowpass(std::array<double, N> const&& data) : filter<Order + 1>
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
            item = (i == half ? 1.0 : 0.0) - factor * filter<Order + 1>::sinc(factor * ((double)(i) - (double)half));
            i++;
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
    explicit highpass(std::array<double, N> const& data) : filter<Order + 1>
    {
        prepare_table(), data
    }
    {}

    /**
     * Constructs a highpass filter with the given cutoff frequency and sampling rate, and applies it to the given data
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit highpass(std::array<double, N> const&& data) : filter<Order + 1>
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
        constexpr std::array<double, Order + 1> low  = lowpass<Order, LowerCutoff, Sampling>::prepare_table();
        constexpr std::array<double, Order + 1> high = highpass<Order, HigherCutoff, Sampling>::prepare_table();
        std::array<double, Order + 1>           array{};
        auto                                    low_ptr = low.begin();
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
    explicit bandstop(std::array<double, N> const& data) : filter<Order + 1>
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
    explicit bandstop(std::array<double, N> const&& data) : filter<Order + 1>
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
    explicit bandpass(std::array<double, N> const& data) : filter<Order + 1>
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
    explicit bandpass(std::array<double, N> const&& data) : filter<Order + 1>
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
    explicit moving_average(std::array<double, N> const& data) : filter<Order>
    {
        prepare_table(), data
    }
    {}

    /**
     * Constructs a moving average filter with the given order and applies it to the given data
     * @param data the data to filter
     */
    template <std::size_t N>
    explicit moving_average(std::array<double, N> const&& data) : filter<Order>
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
