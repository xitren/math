#include <xitren/math/pid2.hpp>

#include <gtest/gtest.h>

#include <iostream>

using pid = xitren::math::pid2_f;

bool
compare_float(float const left, float const right)
{
    constexpr float precision{0.001f};

    if (std::abs(left - right) > precision) {
        std::cout << left << " != " << right << std::endl;
        return false;
    } else {
        return true;
    }
}

class dummy_integrator {
public:
    explicit dummy_integrator(float const ts) : ts_(ts) {}
    dummy_integrator() = delete;

    float
    update(float value)
    {
        auto ret = value_;
        value_ += value * ts_;
        return ret;
    }

    float
    get() const
    {
        return value_;
    }

private:
    float ts_{}, value_{};
};

TEST(pid2_test, pid2_output_proportional)
{
    pid test{0.1f, 0.f, 0.5f, 0.f, 0.f, 10.f, -10.f};

    float target{1.0f};
    for (int i{}; i < 10; i++) {
        test.value(target - 0.f);
    }
    EXPECT_TRUE(compare_float(test.value(target - 0.f), 0.5f));

    target = 5.f;
    for (int i{}; i < 10; i++) {
        test.value(target - 0.f);
    }
    EXPECT_TRUE(compare_float(test.value(target - 0.f), 5 * 0.5f));
}

TEST(pid2_test, pid2_output_integral)
{
    constexpr float ts{0.1f}, ki{0.5f};

    float setpoint{0.5f};
    float control_value{}, next_control_value{};

    pid test{ts, 0.f, 0.f, ki, 0.f, 10.f, -10.f};

    for (int i{}; i < 10; i++) {
        // Backward euler method
        control_value = next_control_value;
        EXPECT_TRUE(compare_float(test.value(setpoint), control_value));
        next_control_value += setpoint * ts * ki;
    }

    // setpoint jump
    setpoint = 2.f;
    for (int i{}; i < 5; i++) {
        control_value = next_control_value;
        EXPECT_TRUE(compare_float(test.value(setpoint), control_value));
        next_control_value += setpoint * ts * ki;
    }
}

TEST(pid2_test, pid2_output_derivative_no_filter)
{
    constexpr float ts{0.1f}, fc_hz{0.f}, kd{1.f};
    float           setpoint{1.f};

    pid test{ts, fc_hz, 0.f, 0.f, kd, 10.f, -10.f};

    float control_value = setpoint * kd / ts;
    EXPECT_TRUE(compare_float(test.value(setpoint), control_value));
}

TEST(pid2_test, pid2_output_proportional_saturation)
{
    constexpr float limit{1.f};

    pid test{0.1f, 0.f, 1.f, 0.f, 0.f, limit, -limit};

    // positive setpoint
    auto target{2.f};
    EXPECT_TRUE(compare_float(test.value(target - (-4.f)), limit));
    EXPECT_TRUE(compare_float(test.value(target - (-2.f)), limit));
    EXPECT_TRUE(compare_float(test.value(target - 0.f), limit));
    EXPECT_TRUE(compare_float(test.value(target - 2.f), 0.f));
    EXPECT_TRUE(compare_float(test.value(target - 4.f), -limit));

    // zero setpoint
    target = 0.f;
    EXPECT_TRUE(compare_float(test.value(target - (-4.f)), limit));
    EXPECT_TRUE(compare_float(test.value(target - (-2.f)), limit));
    EXPECT_TRUE(compare_float(test.value(target - 0.f), 0.f));
    EXPECT_TRUE(compare_float(test.value(target - 2.f), -limit));
    EXPECT_TRUE(compare_float(test.value(target - 4.f), -limit));

    // negative setpoint
    target = -2.f;
    EXPECT_TRUE(compare_float(test.value(target - (-4.f)), limit));
    EXPECT_TRUE(compare_float(test.value(target - (-2.f)), 0.f));
    EXPECT_TRUE(compare_float(test.value(target - 0.f), -limit));
    EXPECT_TRUE(compare_float(test.value(target - 2.f), -limit));
    EXPECT_TRUE(compare_float(test.value(target - 4.f), -limit));
}

TEST(pid2_test, pid2_output_integral_saturation_positive_sp)
{
    constexpr float limit{1.f};

    pid test{0.1f, 0.f, 0.f, 10.f, 0.f, limit, -limit};

    // positive setpoint
    auto target{2.f};

    float y = -2.f;    // positive error
    for (int i{}; i < 10; i++) {
        test.value(target - y);
    }
    EXPECT_TRUE(compare_float(test.value(target - y), limit));

    y = 2.f;    // no error
    for (int i{}; i < 10; i++) {
        test.value(target - y);
    }
    // then get same limited value
    EXPECT_TRUE(compare_float(test.value(target - y), limit));

    y = 4.f;    // negative error
    for (int i{}; i < 10; i++) {
        test.value(target - y);
    }
    EXPECT_TRUE(compare_float(test.value(target - y), -limit));
}

TEST(pid2_test, pid2_output_integral_saturation_negative_sp)
{
    constexpr float limit{1.f};

    pid test{0.1f, 0.f, 0.f, 10.f, 0.f, limit, -limit};

    // negative setpoint
    auto target{-2.f};

    float y = 2.f;    // negative error
    for (int i{}; i < 10; i++) {
        test.value(target - y);
    }
    EXPECT_TRUE(compare_float(test.value(target - y), -limit));

    y = -2.f;    // no error
    for (int i{}; i < 10; i++) {
        test.value(target - y);
    }
    // then get same limited value
    EXPECT_TRUE(compare_float(test.value(target - y), -limit));

    y = -4.f;    // positive error
    for (int i{}; i < 10; i++) {
        test.value(target - y);
    }
    EXPECT_TRUE(compare_float(test.value(target - y), limit));
}

TEST(pid2_test, pid2_integral_plant_proportional_output)
{
    pid              test{0.1f, 0.f, 1.f, 0.f, 0.f, 10.f, -10.f};
    dummy_integrator plant{0.1f};

    auto target{1.f};
    for (int i{}; i < 100; i++) {
        plant.update(test.value(target - plant.get()));
    }
    EXPECT_TRUE(compare_float(plant.get(), 1.f));
}

TEST(pid2_test, pid2_integral_plant_proportional_output_saturated)
{
    pid              test{0.1f, 0.f, 1.f, 0.f, 0.f, 0.5f, -0.5f};
    dummy_integrator plant{0.1f};

    auto target{1.f};
    for (int i{}; i < 100; i++) {
        plant.update(test.value(target - plant.get()));
    }
    EXPECT_TRUE(compare_float(plant.get(), 1.f));
}

// Expects rising of plant value
TEST(pid2_test, pid2_integral_plant_integral_output)
{
    pid              test{0.1f, 0.f, 0.f, 0.1f, 0.f, 10.f, -10.f};
    dummy_integrator plant{0.1f};

    auto target{1.f};

    bool ok{false};
    for (int i{}; i < 200; i++) {
        plant.update(test.value(target - plant.get()));
        if (plant.get() >= 1.f) {
            ok = true;
            break;
        }
    }
    EXPECT_TRUE(ok);
}

TEST(pid2_test, pid2_integral_plant_integral_output_saturated)
{
    pid              test{0.1f, 0.f, 0.f, 0.1f, 0.f, 0.5f, 0.5f};
    dummy_integrator plant{0.1f};

    auto target{1.f};

    bool ok{false};
    for (int i{}; i < 200; i++) {
        plant.update(test.value(target - plant.get()));
        if (plant.get() >= 1.f) {
            ok = true;
            break;
        }
    }
    EXPECT_TRUE(ok);
}

TEST(pid2_test, pid2_integral_reset)
{
    pid test{0.1f, 0.f, 0.f, 1.f, 0.f, 10.f, -10.f};

    auto target{1.f};

    // accumulate integral part
    test.value(target);
    test.value(target);
    EXPECT_FALSE(compare_float(test.value(target), 0));

    // reset integral part
    // will be zero only for backward-euler
    test.reset();
    EXPECT_TRUE(compare_float(test.value(target), 0));
}
