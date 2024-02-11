#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <numbers>

namespace xitren::math {

template <typename T>
concept rational_type = std::same_as<T, double> || std::same_as<T, float>;

template <rational_type T = double>
struct pid_params {
    T sampling_time{1.f};
    T filter{20.f};
    T kp{1.f};
    T ki{1.f};
    T kd{1.f};
    T max{1000.f};
    T min{-1000.f};
};

template <rational_type rational_t>
class pid2 {

    // Compute the weight factor alpha for an exponential moving average LPF filter
    // with a given normalized cutoff frequency `fn` = f_cutoff_Hz  * t_sampling_s.
    static constexpr rational_t
    update_alpha(const rational_t fn)
    {
        constexpr rational_t min_freq_hz{1.f};
        constexpr auto       pi = std::numbers::pi_v<rational_t>;

        if (fn < min_freq_hz)
            return 1;

        // alpha(frq) = cos(2 * Pi * frq) - 1
        // + sqrt{ cos(2 * Pi * frq)^2 - 4 * cos(2 * Pi * frq) + 3}

        // let c = cos(2 * Pi * frq), then:
        // alpha = c - 1 + sqrt( c * c - 4 * c + 3)

        const rational_t c = std::cos(2 * pi * fn);
        return c - 1 + std::sqrt(c * c - 4 * c + 3);
    };

public:
    // ts  Controller sampling time in seconds
    // fc  Cutoff frequency of derivative filter in Hz
    // kp  Proportional gain
    // ki  Integral gain
    // kd  Derivative gain

    // The derivative filter can be disabled by setting `fc` less than 1 Hz.

    constexpr explicit pid2(rational_t ts, rational_t fc_hz, rational_t kp, rational_t ki, rational_t kd,
                            rational_t max, rational_t min)
        : ts_(ts), fc_hz_(fc_hz), kp_(kp), ki_(ki), kd_(kd), alpha_(update_alpha(fc_hz_ * ts_)), max_(max), min_(min)

    {}

    constexpr explicit pid2(pid_params<rational_t> const& p)
        : pid2(p.sampling_time, p.filter, p.kp, p.ki, p.kd, p.max, p.min)
    {}

    // Get the controller output with the given 'error'
    rational_t
    value(rational_t error)
    {
        // Input is e[k] = r[k] - y[k], error between setpoint and true position

        // e_f[k] = alpha * e[k] + (1-alpha) * e_f[k-1], filtered error
        rational_t ef = alpha_ * error + (1 - alpha_) * ef_prev_;

        // e_d[k] = (e_f[k] - e_f[k-1]) / Ts, filtered derivative part
        rational_t derivative = (ef - ef_prev_) / ts_;

        // e_i[k+1] = e_i[k] + Ts * e[k], integral part for next time
        rational_t next_integral_value = integral_value_ + error * ts_;

        // PID formula:
        // u[k] = Kp * e[k] + Ki * e_i[k] + Kd * e_d[k]
        rational_t control_u = kp_ * error + ki_ * integral_value_ + kd_ * derivative;

        // Limit output
        if (control_u > max_) {
            control_u = max_;
            // check sign of integral value delta to escape from saturation
            if (next_integral_value < integral_value_)
                integral_value_ = next_integral_value;
        } else if (control_u < min_) {
            control_u = min_;
            // check sign of integral value delta to escape from saturation
            if (next_integral_value > integral_value_)
                integral_value_ = next_integral_value;
        } else {    // anti-windup
            integral_value_ = next_integral_value;
        }

        // store the state for the next iteration
        ef_prev_ = ef;
        return control_u;
    }

    void
    reset()
    {
        ef_prev_        = 0;
        integral_value_ = 0;
    }

    [[nodiscard]] rational_t
    max() const
    {
        return max_;
    }

    pid2&
    max(rational_t val)
    {
        max_ = val;
        return *this;
    }

    [[nodiscard]] rational_t
    min() const
    {
        return min_;
    }

    pid2&
    min(rational_t val)
    {
        min_ = val;
        return *this;
    }

    [[nodiscard]] rational_t
    proportional() const
    {
        return kp_;
    }

    pid2&
    proportional(rational_t val)
    {
        kp_ = val;
        return *this;
    }

    [[nodiscard]] rational_t
    integral() const
    {
        return ki_;
    }

    pid2&
    integral(rational_t val)
    {
        ki_ = val;
        return *this;
    }

    rational_t
    get_integrated()
    {
        return integral_value_;
    }

    [[nodiscard]] rational_t
    derivative() const
    {
        return kd_;
    }

    pid2&
    derivative(rational_t val)
    {
        kd_ = val;
        return *this;
    }

    [[nodiscard]] rational_t
    filter() const
    {
        return fc_hz_;
    }

    pid2&
    filter(rational_t val)
    {
        fc_hz_ = val;
        update_alpha(fc_hz_ * ts_);
        return *this;
    }

    [[nodiscard]] rational_t
    sampling_time() const
    {
        return ts_;
    }
    pid2&
    sampling_time(rational_t val)
    {
        ts_ = val;
        update_alpha(fc_hz_ * ts_);
        return *this;
    }

private:
    rational_t ts_, fc_hz_, kp_, ki_, kd_, alpha_, max_, min_;
    rational_t integral_value_ = 0;
    rational_t ef_prev_        = 0;
};

using pid2_f = pid2<float>;
using pid2_d = pid2<double>;

using pid_params_f = pid_params<float>;
using pid_params_d = pid_params<double>;

}    // namespace loveka::components::math
