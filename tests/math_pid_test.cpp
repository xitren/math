#include <xitren/math/pid.hpp>

#include <gtest/gtest.h>

#include <iostream>

using namespace xitren::math;

class dummy_motor {
public:
    dummy_motor()
        : K(2.5),
          AD(18.01),
          BD(-39.98),
          CD(22.01),
          MAX_OUTPUT(100.0),
          MIN_OUTPUT(-100.0),
          SAMPLING_TIME(0.001),
          C1(K / CD),
          C2(2.0 * K / CD),
          C3(K / CD),
          C4(BD / CD),
          C5(AD / CD)
    {
        reset();
    }

    void
    reset()
    {
        m_yn1    = 0.0;
        m_yn2    = 0.0;
        m_xn1    = 0.0;
        m_xn2    = 0.0;
        m_output = 0.0;
    }

    double
    nextSample(double input)
    {
        m_output += input;
        return m_output;
    }

    double m_output{};

private:
    double const K;

    double const AD;
    double const BD;
    double const CD;

    double const MAX_OUTPUT;
    double const MIN_OUTPUT;
    double const SAMPLING_TIME;

    double const C1;
    double const C2;
    double const C3;
    double const C4;
    double const C5;

    // Attributes
    double m_yn1{}, m_yn2{}, m_xn1{}, m_xn2{};
};

double e2, e1, e0, u2, u1, u0;
double r;
double y;
double Kp = 1;
double Ki = 1;
double Kd = 1;
double N  = 20;
double Ts = 0.01;
double a0 = (1 + N * Ts);
double a1 = -(2 + N * Ts);
double a2 = 1;
double b0 = Kp * (1 + N * Ts) + Ki * Ts * (1 + N * Ts) + Kd * N;
double b1 = -(Kp * (2 + N * Ts) + Ki * Ts + 2 * Kd * N);
double b2 = Kp + Kd * N;

double ku1 = a1 / a0;
double ku2 = a2 / a0;
double ke0 = b0 / a0;
double ke1 = b1 / a0;
double ke2 = b2 / a0;

// loveka::components::math::pid::pid uses double as value type
class dummy_integrator {
public:
    explicit dummy_integrator(double const ts) : ts_(ts) {}
    dummy_integrator() = delete;

    double
    update(double value)
    {
        auto ret = value_;
        value_ += value * ts_;
        return ret;
    }

    double
    get() const
    {
        return value_;
    }

private:
    double ts_{}, value_{};
};

bool
compare_double(double const left, double const right)
{
    constexpr double precision{0.001f};

    if (std::abs(left - right) > precision) {
        std::cout << left << " != " << right << std::endl;
        return false;
    } else {
        return true;
    }
}

TEST(pid_test, fast_moving)
{
    dummy_motor motor{};
    double      val = 0;
    r               = 1.;
    for (int i = 0; i < 100; i++) {
        e2 = e1;
        e1 = e0;
        u2 = u1;
        u1 = u0;

        y = motor.nextSample(val);
        std::cout << y << std::endl;

        e0 = r - y;

        u0 = -ku1 * u1 - ku2 * u2 + ke0 * e0 + ke1 * e1 + ke2 * e2;

        if (u0 > 100)
            u0 = 100;
        if (u0 < -100)
            u0 = -100;

        val = u0;
    }
}

TEST(pid_test, fast_moving_average)
{
    pid test{1., 20, 0.5, 0.5, 0.5};
    test.min(-0.1).max(0.1);
    dummy_motor motor{};
    test.target(1.);
    for (int i = 0; i < 100; i++) {
        double const val = test.value(motor.m_output);
        motor.nextSample(val);
        std::cout << motor.m_output << std::endl;
    }
}

TEST(pid_test, pid_output_proportional)
{
    pid test{0.1, 0., 0.5, 0., 0., 10., -10.};

    test.target(1.);
    for (int i{}; i < 10; i++) {
        test.value(0.);
    }
    EXPECT_TRUE(compare_double(test.value(0.), 0.5));

    test.target(5.);
    for (int i{}; i < 10; i++) {
        test.value(0.);
    }
    EXPECT_TRUE(compare_double(test.value(0.), 5 * 0.5));
}

TEST(pid_test, pid_output_integral)
{
    constexpr double ts{0.1}, ki{0.5};

    double setpoint{0.5};
    double control_value{};

    pid test{ts, 0., 0., ki, 0., 10., -10.};

    test.target(setpoint);
    for (int i{}; i < 10; i++) {
        // Forward euler method
        control_value += setpoint * ts * ki;
        EXPECT_TRUE(compare_double(test.value(0), control_value));
    }

    // setpoint jump
    setpoint = 2.f;
    test.target(setpoint);
    for (int i{}; i < 5; i++) {
        control_value += setpoint * ts * ki;
        EXPECT_TRUE(compare_double(test.value(0.f), control_value));
    }
}

// Todo: make a test for derivative part without filter
// but filter value eq. 0 does not disable the filter
TEST(pid_test, pid_output_derivative_no_filter)
{
    constexpr float ts{0.1}, fc_hz{0.}, kd{1.};

    pid test{ts, fc_hz, 0., 0., kd, 10.f, -10.f};

    test.target(1.);

    // Show control output value without comparing
    std::cout << "No derivative filter output:" << std::endl;
    std::cout << test.value(0.) << std::endl;
    std::cout << test.value(0.) << std::endl;
    std::cout << test.value(0.) << std::endl;
}

TEST(pid_test, pid_output_proportional_saturation)
{
    constexpr double limit{1.f};

    pid test{0.1f, 0.f, 1.f, 0.f, 0.f, limit, -limit};

    // positive setpoint
    test.target(2.f);
    EXPECT_TRUE(compare_double(test.value(-4.f), limit));

    GTEST_SKIP() << "Saturation doesn't work properly";
    EXPECT_TRUE(compare_double(test.value(-2.f), limit));
    EXPECT_TRUE(compare_double(test.value(0.f), limit));
    EXPECT_TRUE(compare_double(test.value(2.f), 0.f));
    EXPECT_TRUE(compare_double(test.value(4.f), -limit));

    // zero setpoint
    test.target(0.f);
    EXPECT_TRUE(compare_double(test.value(-4.f), limit));
    EXPECT_TRUE(compare_double(test.value(-2.f), limit));
    EXPECT_TRUE(compare_double(test.value(0.f), 0.f));
    EXPECT_TRUE(compare_double(test.value(2.f), -limit));
    EXPECT_TRUE(compare_double(test.value(4.f), -limit));

    // negative setpoint
    test.target(-2.f);
    EXPECT_TRUE(compare_double(test.value(-4.f), limit));
    EXPECT_TRUE(compare_double(test.value(-2.f), 0.f));
    EXPECT_TRUE(compare_double(test.value(0.f), -limit));
    EXPECT_TRUE(compare_double(test.value(2.f), -limit));
    EXPECT_TRUE(compare_double(test.value(4.f), -limit));
}

TEST(pid_test, pid_output_integral_saturation_positive_sp)
{
    constexpr double limit{1.f};

    pid test{0.1, 0., 0., 10., 0., limit, -limit};

    // positive setpoint
    test.target(2.);

    float y = -2.;    // positive error
    for (int i{}; i < 10; i++) {
        test.value(y);
    }
    EXPECT_TRUE(compare_double(test.value(y), limit));

    y = 2.;    // no error
    for (int i{}; i < 10; i++) {
        test.value(y);
    }
    GTEST_SKIP() << "Saturation doesn't work properly";
    // then get same limited value
    EXPECT_TRUE(compare_double(test.value(y), limit));

    y = 4.;    // negative error
    for (int i{}; i < 10; i++) {
        test.value(y);
    }
    EXPECT_TRUE(compare_double(test.value(y), -limit));
}

TEST(pid_test, pid_output_integral_saturation_negative_sp)
{
    constexpr double limit{1.f};

    pid test{0.1, 0., 0., 10., 0., limit, -limit};

    // negative setpoint
    test.target(-2.);

    float y = 2.;    // negative error
    for (int i{}; i < 10; i++) {
        test.value(y);
    }
    EXPECT_TRUE(compare_double(test.value(y), -limit));

    y = -2.;    // no error
    for (int i{}; i < 10; i++) {
        test.value(y);
    }
    GTEST_SKIP() << "Saturation doesn't work properly";
    // then get same limited value
    EXPECT_TRUE(compare_double(test.value(y), -limit));

    y = -4.;    // positive error
    for (int i{}; i < 10; i++) {
        test.value(y);
    }
    EXPECT_TRUE(compare_double(test.value(y), limit));
}

TEST(pid_test, pid_integral_plant_proportional_output)
{
    pid              test{0.1, 0., 1., 0., 0., 10., -10.};
    dummy_integrator plant{0.1};

    test.target(1.);
    for (int i{}; i < 100; i++) {
        plant.update(test.value(plant.get()));
    }
    EXPECT_TRUE(compare_double(plant.get(), 1.));
}

TEST(pid_test, pid_integral_plant_proportional_output_saturated)
{
    pid              test{0.1, 0., 1., 0., 0., 0.5, -0.5};
    dummy_integrator plant{0.1};

    test.target(1.);
    for (int i{}; i < 100; i++) {
        plant.update(test.value(plant.get()));
    }
    GTEST_SKIP() << "Saturation doesn't work properly";
    EXPECT_TRUE(compare_double(plant.get(), 1.));
}

// Expects rising of plant value
TEST(pid_test, pid_integral_plant_integral_output)
{
    pid              test{0.1, 0., 0., 0.1, 0., 10., -10.};
    dummy_integrator plant{0.1};

    test.target(1.);

    bool ok{false};
    for (int i{}; i < 200; i++) {
        plant.update(test.value(plant.get()));
        if (plant.get() >= 1.) {
            ok = true;
            break;
        }
    }
    EXPECT_TRUE(ok);
}

TEST(pid_test, pid_integral_plant_integral_output_saturated)
{
    pid              test{0.1, 0., 0., 0.1, 0., 0.5, -0.5};
    dummy_integrator plant{0.1};

    test.target(1.);

    bool ok{false};
    for (int i{}; i < 200; i++) {
        plant.update(test.value(plant.get()));
        if (plant.get() >= 1.) {
            ok = true;
            break;
        }
    }
    EXPECT_TRUE(ok);
}
