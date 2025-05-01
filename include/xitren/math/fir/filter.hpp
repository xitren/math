/*!
_ _
__ _(_) |_ _ _ ___ _ _
\ \ / |  _| '_/ -_) ' \
/_\_\_|\__|_| \___|_||_|
* @date 20.01.2025
*/
#pragma once

#include <xitren/circular_buffer.hpp>

#include <array>
#include <cmath>
#include <cstdint>

namespace xitren::math::fir {

template <typename Type, std::size_t Order, typename TableType = double>
class filter : public circular_buffer<Type, Order> {
    static_assert(std::is_floating_point<TableType>(), "Table should contain floating point values!");

    using circular_buffer<Type, Order>::begin;
    using circular_buffer<Type, Order>::end;
    using circular_buffer<Type, Order>::full;

    static constexpr double M_2PI = 2. * M_PI;

protected:
    using data_type   = std::array<Type, Order>;
    using table_type  = std::array<TableType, Order>;
    using filter_type = filter;

public:
    /**
     * Constructs a filter with the given table data
     * @param table_data the table data to use for the filter
     */
    constexpr explicit filter(table_type const& table_data) : table_{table_data} {}

    /**
     * Constructs a filter with the given table data and data
     * @param table_data the table data to use for the filter
     * @param data the data to filter
     */
    constexpr filter(table_type const& table_data, data_type const& data) : table_{table_data} { (*this) << data; }

    /**
     * Constructs a filter with the given table data and rvalue data
     * @param table_data the table data to use for the filter
     * @param data the rvalue data to filter
     */
    filter(table_type const& table_data, data_type const&& data) : table_{table_data} { (*this) << data; }

    consteval filter
    operator*(filter const& other) const
    {
        filter ret(*this);
        for (std::size_t i{}; i < Order; ++i) {
            ret.table_[i] = this->table_[i] * other.table_[i];
        }
        return ret;
    }

    consteval filter
    operator+(filter const& other) const
    {
        filter ret(*this);
        for (std::size_t i{}; i < Order; ++i) {
            ret.table_[i] = this->table_[i] * other.table_[i];
        }
        return ret;
    }

    consteval filter
    operator-(filter const& other) const
    {
        filter ret(*this);
        for (std::size_t i{}; i < Order; ++i) {
            ret.table_[i] = this->table_[i] * other.table_[i];
        }
        return ret;
    }

    /**
     * Applies the filter to a new data point
     * @param val the new data point
     * @return the filtered data point
     */
    Type
    value(Type val)
    {
        circular_buffer<Type, Order>::push(val);
        if (!full())
            return 0.;
        auto      it = begin();
        TableType ret_val{};
        for (auto& item : table_) {
            ret_val += item * (*it);
            ++it;
        }
        return static_cast<Type>(ret_val);
    }

    /**
     * Resets the filter state
     */
    void
    reset()
    {
        circular_buffer<Type, Order>::clear();
    }

    /**
     * Returns the filter table data
     * @return the filter table data
     */
    std::array<TableType, Order>
    table() const
    {
        return table_;
    }

private:
    table_type table_;

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
    sin_cfrac(Real x2, int const k = 2, int const n = 40)
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
        return (x <= -M_PI) ? wrap(x + M_2PI) : (x > M_PI) ? wrap(x - M_2PI) : (true) ? x : 0;
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
        if (x == 0) {
            return 1.0;
        }
        double const xpi = M_PI * x;
        return sin(xpi) / xpi;
    }

    /**
     * Calculates the modified Bessel function of the first kind of order 0
     * @param x the argument of the function
     * @return the modified Bessel function of the first kind of order 0 evaluated at x
     */
    static constexpr double
    i0(double const x)
    {
        double       f  = 1.;
        double const x2 = x * x * 0.25;
        double       xc = x2;
        double       v  = 1. + x2;
        for (int i = 2; i < 100; i++) {
            f *= static_cast<double>(i);
            xc *= x2;
            double const a = xc / (f * f);
            v += a;
            if (a < 1e-20) {
                break;
            }
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
        double sum{};
        for (auto& item : win) {
            sum += item;
        }
        if (sum == 0.) {
            return;
        }
        for (auto& item : win) {
            item /= sum;
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
        std::int32_t m  = (attenuation <= 21.) ? static_cast<std::int32_t>(::ceil(5.79 / tw))
                                               : static_cast<std::int32_t>(::ceil((attenuation - 7.95) / (2.285 * tw)));
        if ((m & 1) == 0) {
            m++;
        }
        double const beta
            = (attenuation <= 21.)
                  ? (0.)
                  : ((attenuation <= 50.) ? (0.5842 * ::pow(attenuation - 21., 0.4) + 0.07886 * (attenuation - 21.))
                                          : (0.1102 * (attenuation - 8.7)));
        double const i0b = i0(beta);
        for (int n = 0; n < m; n++) {
            double const v
                = beta * ::sqrt(1.0 - ::pow(2.0 * static_cast<double>(n) / (static_cast<double>(m) - 1.) - 1.0, 2.));
            win[n] = i0(v) / i0b;
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
        double const tw = 2.0 * M_PI / (static_cast<double>(Size) - 1.);
        for (std::size_t i{}; i < Size; i++) {
            win[i] *= 0.42 - 0.5 * ::cos(tw * static_cast<double>(i)) + 0.08 * ::cos(2.0 * tw * static_cast<double>(i));
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
        std::size_t       i{};
        for (auto& item : win) {
            item *= sinc(2.0 * static_cast<double>(i++) / static_cast<double>(m) - 1.0);
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
        std::size_t       i{};
        for (auto& item : win) {
            item *= 0.5 - 0.5 * ::cos(2.0 * M_PI * static_cast<double>(i++) / static_cast<double>(m));
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
        std::size_t       i{};
        for (auto& item : win) {
            item *= 0.54 - 0.46 * ::cos(2.0 * M_PI * static_cast<double>(i++) / static_cast<double>(m));
        }
    }
};

}    // namespace xitren::math::fir
