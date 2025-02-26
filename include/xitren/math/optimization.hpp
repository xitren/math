#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <functional>
#include <memory>
#include <thread>
#include <utility>
#include <vector>
using namespace std::literals;

namespace xitren::math {

template <class Type, std::size_t Inputs, std::size_t Outputs>
class optimization {
    using time_type     = std::chrono::milliseconds;
    using input_type    = std::array<Type, Inputs>;
    using output_type   = std::array<Type, Outputs>;
    using run_func_type = std::function<output_type(input_type)>;

public:
    optimization() = delete;
    optimization&
    operator=(optimization const& other)
        = delete;
    optimization&
    operator=(optimization const&& other)
        = delete;
    optimization(optimization const& val) = delete;
    optimization(optimization&& val)      = delete;

    static double
    quadratic(output_type const& target, output_type const& current)
    {
        double sum{};
        for (std::size_t i{}; i < target.size(); i++) {
            sum += std::abs(current[i] - target[i]);
        }
        return sum / target.size();
    }

    input_type
    gradient(input_type const& current, double lamda) const
    {
        input_type grad{};
        input_type dis{current};
        for (std::size_t i{}; i < current.size(); i++) {
            dis[i] += lamda;
            auto ip   = function_(dis);
            auto errp = quadratic(target_, ip);
            dis[i] -= lamda;
            dis[i] -= lamda;
            auto in   = function_(dis);
            auto errn = quadratic(target_, in);
            dis[i] += lamda;
            grad[i] = (errp - errn) / (2 * lamda);
        }
        return grad;
    }

    double
    get_lambda(input_type const& grad, double step)
    {
        double lambda{step};
        auto   current{current_};
        for (std::size_t i{}; i < current_.size(); i++) {
            current[i] -= lambda * grad[i];
        }
        double old{quadratic(target_, function_(current))};
        double new_one{old - 1};
        for (; old > new_one;) {
            lambda *= 2;
            auto next{current_};
            for (std::size_t i{}; i < current_.size(); i++) {
                next[i] -= lambda * grad[i];
            }
            new_one = quadratic(target_, function_(next));
        }
        return lambda / 2;
    }

    /**
     * Constructs a new optimization object.
     *
     * @param function The function to be called periodically.
     * @param maximum_time The interval at which the function should be called, in milliseconds.
     * @param target target data.
     * @param precision precision.
     */
    optimization(run_func_type const& function, output_type const& target, input_type const& initial, double precision,
                 time_type maximum_time = 100ms)
        : function_{std::move(function)},
          target_{target},
          precision_{std::abs(precision)},
          current_{initial},
          thread_{[&] {
              auto       last_time{std::chrono::high_resolution_clock::now()};
              double     err{};
              double     lambda = 0.01;
              input_type grad{};
              while (keep_running_.test_and_set()) {
                  auto ret     = function_(current_);
                  auto old_err = err;
                  err          = quadratic(target_, ret);
                  auto pass    = std::chrono::duration_cast<std::chrono::milliseconds>(
                                  std::chrono::high_resolution_clock::now() - last_time)
                                  .count();
                  auto expect = (maximum_time_).count();
                  if (pass >= expect || (err <= precision_) || (lambda != lambda)
                      || (std::abs(old_err - err) <= lambda)) {
                      keep_running_.clear();
                      break;
                  }
                  // Calculate gradient
                  auto grad = gradient(current_, lambda);
                  // Calculate lambda
                  lambda = get_lambda(grad, 0.001);
                  // Calculate next point
                  for (std::size_t i{}; i < current_.size(); i++) {
                      current_[i] -= lambda * grad[i];
                  }
              }
              keep_running_.notify_one();
          }},
          maximum_time_{maximum_time}
    {}

    /**
     * Destroys the optimization object.
     */
    ~optimization() { stop(); }

    /**
     * Stops the interval_event object.
     */
    void
    stop()
    {
        if (keep_running_.test_and_set()) {
            // keep_running_.clear();
            keep_running_.wait(false);
        }
        if (thread_.joinable()) {
            thread_.join();
        }
    }

    /**
     * Returns a reference to the thread object that is used to run the interval_event.
     */
    auto&
    thread() noexcept
    {
        return thread_;
    }

    /**
     * Returns a reference to the thread object that is used to run the interval_event.
     */
    auto&
    keep_running() const noexcept
    {
        return keep_running_;
    }

    /**
     * Returns the maximum time at which the function is called, in milliseconds.
     */
    auto&
    maximum_time() const noexcept
    {
        return maximum_time_;
    }

    /**
     * Sets the maximum time at which the function is called.
     *
     * @param val The new interval, in milliseconds.
     */
    void
    maximum_time(time_type const& val) noexcept
    {
        maximum_time_ = val;
    }

    auto&
    result() const noexcept
    {
        return current_;
    }

private:
    run_func_type    function_;
    std::atomic_flag keep_running_{true};
    std::thread      thread_;
    time_type        maximum_time_{};
    output_type      target_;
    input_type       current_;
    double           precision_;
};
}    // namespace xitren::math
