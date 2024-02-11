#pragma once
#include <xitren/math/pid2.hpp>

#include <algorithm>
#include <cstdint>

namespace xitren::math {

using pid_params2 = xitren::math::pid_params<double>;

/*
 * Related article: https://nda.ya.ru/t/wHPMbLAn6L25zE
 * PID:
 *      u(t) = Kp*e(t) + Ki*|e(t)dt + Kd(de(t)/dt)
 *Laplace transform PID:
 *      C(s) = Kp + Ki / s + Kd * s
 *Laplace transform PID with derivative filter:
 *      C(s) = Kp + Ki / s + (N * Kd) / (1 + N * s)
 *Backward Euler integral:
 *      Ki / s = Ki * Ts / (1 - z[-1])
 *Backward Euler derivative:
 *      (N * Kd) / (1 + N * s) = N * (1 - z[-1]) / (1 + N * Ts) - z[-1])
 *Z-domain PID:
 *      U(z) = C(z) * E(z)
 *      C(z) = Kp + Ki * Ts / (1 - z[-1]) + Kd * (N * (1 - z[-1]) / (1 + N * Ts - z[-1])))
 *Z-domain discrete PID:
 *      C(z) = U(z)/E(z) = (b0 + b1 * z[-1] + b2 * z[-2])/(a0 + a1 * z[-1] + a2 * z[-2])
 *          b0 = Kp * (1 + N * Ts) + Ki * Ts * (1 + N * Ts) + Kd * N
 *          b1 = -(Kp * (2 + N * Ts) + Ki * Ts + 2 * Kd * N)
 *          b2 = Kp + Kd * N
 *          a0 = (1 + N * Ts)
 *          a1 = -(2 + N * Ts)
 *          a2 = 1
 *
 * Rearrange:
 *      a0 * U(z) + a1 * z[-1] * U(z) + a2 * z[-2] * U(z) = b0 * E(z) + b1 * z[-1] * E(z) + b2 *
 *z[-2] * E(z) a0 * U(z) = -a1 * z[-1] * U(z) - a2 * z[-2] * U(z) + b0 * E(z) + b1 * z[-1] * E(z) +
 *b2 * z[-2] * E(z) u[k] = -(a1/a0) * u[k - 1] - (a2/a0) * u[k - 2] + (b0/a0) * e[k] + (b1/a0) * e[k
 *- 1] + (b2/a0) * u[k - 2]
 *
 *  Then controller output:
 *      u[k] = -(a1/a0) * u[k - 1] - (a2/a0) * u[k - 2] + (b0/a0) * e[k] + (b1/a0) * e[k - 1] +
 *  + (b2/a0) * e[k - 2]
 *
 * For easily programming describe new coefficients:
 *  ku1 = a1/a0
 *  ku2 = a2/a0
 *  ke0 = b0/a0
 *  ke1 = b1/a0
 *  ke2 = b2/a0
 **/

class pid {
public:
    constexpr explicit pid(pid_params2 const& p) : pid(p.sampling_time, p.filter, p.kp, p.ki, p.kd, p.max, p.min) {}
    constexpr explicit pid(double sampling_time = 1., double filter = 20., double kp = 20., double ki = 1.,
                           double kd = 1., double max = 10000., double min = -10000.)
        : max_(max),
          min_(min),
          proportional_{kp},
          integral_{ki},
          derivative_{kd},
          filter_{filter},
          sampling_time_{sampling_time},
          ku1_{ku1(filter, sampling_time)},
          ku2_{ku2(filter, sampling_time)},
          ke0_{ke0(kp, ki, kd, filter, sampling_time)},
          ke1_{ke1(kp, ki, kd, filter, sampling_time)},
          ke2_{ke2(kp, ki, kd, filter, sampling_time)}
    {}

    double
    value(double y)
    {
        e2_ = e1_;
        e1_ = e0_;
        u2_ = u1_;
        u1_ = u0_;
        e0_ = target_ - y;
        u0_ = (-ku1_ * u1_) + (-ku2_ * u2_) + ke0_ * e0_ + ke1_ * e1_ + ke2_ * e2_;
        u0_ = std::min(u0_, max_);
        u0_ = std::max(u0_, min_);
        return u0_;
    }

    // experiment
    void
    reset()
    {
        e2_ = 0;
        e1_ = 0;
        e0_ = 0;
        u2_ = 0;
        u1_ = 0;
        u0_ = 0;
    }

    [[nodiscard]] double
    target() const
    {
        return target_;
    }
    pid&
    target(double val)
    {
        target_ = val;
        return *this;
    }

    [[nodiscard]] double
    max() const
    {
        return max_;
    }
    pid&
    max(double val)
    {
        max_ = val;
        return *this;
    }

    [[nodiscard]] double
    min() const
    {
        return min_;
    }
    pid&
    min(double val)
    {
        min_ = val;
        return *this;
    }

    [[nodiscard]] double
    proportional() const
    {
        return proportional_;
    }
    pid&
    proportional(double val)
    {
        proportional_ = val;
        recalculate();
        return *this;
    }

    [[nodiscard]] double
    integral() const
    {
        return integral_;
    }
    pid&
    integral(double val)
    {
        integral_ = val;
        recalculate();
        return *this;
    }

    [[nodiscard]] double
    derivative() const
    {
        return derivative_;
    }
    pid&
    derivative(double val)
    {
        derivative_ = val;
        recalculate();
        return *this;
    }

    [[nodiscard]] double
    filter() const
    {
        return filter_;
    }
    pid&
    filter(double val)
    {
        filter_ = val;
        recalculate();
        return *this;
    }

    [[nodiscard]] double
    sampling_time() const
    {
        return sampling_time_;
    }
    pid&
    sampling_time(double val)
    {
        sampling_time_ = val;
        recalculate();
        return *this;
    }

private:
    double target_{};
    double max_{10000.};
    double min_{-10000.};
    double proportional_{20.};
    double integral_{1.};
    double derivative_{1.};
    double filter_{20.};
    double sampling_time_ = 0.01;
    double ku1_, ku2_, ke0_, ke1_, ke2_;
    double e2_{}, e1_{}, e0_{}, u2_{}, u1_{}, u0_{};

    void
    recalculate()
    {
        ku1_ = ku1(filter_, sampling_time_);
        ku2_ = ku2(filter_, sampling_time_);
        ke0_ = ke0(proportional_, integral_, derivative_, filter_, sampling_time_);
        ke1_ = ke1(proportional_, integral_, derivative_, filter_, sampling_time_);
        ke2_ = ke2(proportional_, integral_, derivative_, filter_, sampling_time_);
    }

    static constexpr double
    a0(double filter, double sampling_time)
    {
        return 1 + filter * sampling_time;
    }

    static constexpr double
    a1(double filter, double sampling_time)
    {
        return -(2 + filter * sampling_time);
    }

    static constexpr double
    a2()
    {
        return 1;
    }

    static constexpr double
    b0(double kp, double ki, double kd, double filter, double sampling_time)
    {
        return kp * (1 + filter * sampling_time) + ki * sampling_time * (1 + filter * sampling_time) + kd * filter;
    }

    static constexpr double
    b1(double kp, double ki, double kd, double filter, double sampling_time)
    {
        return -(kp * (2 + filter * sampling_time) + ki * sampling_time + 2 * kd * filter);
    }

    static constexpr double
    b2(double kp, double kd, double filter)
    {
        return kp + kd * filter;
    }

    static constexpr double
    ku1(double filter, double sampling_time)
    {
        return a1(filter, sampling_time) / a0(filter, sampling_time);
    }

    static constexpr double
    ku2(double filter, double sampling_time)
    {
        return a2() / a0(filter, sampling_time);
    }

    static constexpr double
    ke0(double kp, double ki, double kd, double filter, double sampling_time)
    {
        return b0(kp, ki, kd, filter, sampling_time) / a0(filter, sampling_time);
    }

    static constexpr double
    ke1(double kp, double ki, double kd, double filter, double sampling_time)
    {
        return b1(kp, ki, kd, filter, sampling_time) / a0(filter, sampling_time);
    }

    static constexpr double
    ke2(double kp, [[maybe_unused]] double ki, double kd, double filter, double sampling_time)
    {
        return b2(kp, kd, filter) / a0(filter, sampling_time);
    }
};

}    // namespace loveka::components::math::pid
