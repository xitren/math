#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <utility>

/** De Casteljau Algoritm
 *
 * Main concept is calculate a lines between base points
 * to get new lines with lower degree.
 * Calculate it repeatedly until last point will be rest
 *
 * 0 <= t <= 1 -- position on curve
 *
 * 3-degree example
 * P0_0(x,y)...P3_0(x,y) - 4 base point
 *
 * Get 3 point next degree:
 * P0_1 = P0_0 + (P1_0 - P0_0) * t
 * P1_1 = P1_0 + (P2_0 - P1_0) * t
 * P2_1 = P2_0 + (P3_0 - P1_0) * t
 *
 * Next step has only 2 point:
 * P0_2 = P0_1 + (P1_1 - P0_1) * t
 * P1_2 = P1_1 + (P2_1 - P0_1) * t
 *
 * And last step to get point on curve:
 * P0_3 = P0_2 + (P1_2 - P0_2) * t
 *
 * If put values of each points to last equality, we get
 * qubic formula:
 * P(t) = k3 * t^3 + k2 * t^2 + k1 * t + k0
 * k3 = -P0 + 3 * P1 - 3 * P2 + P3
 * k2 = 3 * P0 - 6 * P1 + 3 * P2
 * k1 = -P0 + 3 * P1
 * k0 = P0
 *
 */

namespace xitren::math {

enum class degree { quadratic = 3, cubic = 4 };

template <typename T>
concept coordinate_typename
    = std::same_as<float, T> || std::same_as<double, T> || std::same_as<std::int32_t, T>;

template <degree degree_, coordinate_typename T>
class bezier {
public:
    static_assert(degree_ == degree::cubic, "Only cubic bezier is supported yet");

    using coord_t = T;

    struct point_t {
        coord_t x{};
        coord_t y{};
    };

    using points_array_t = std::array<point_t, static_cast<std::uint8_t>(degree_)>;

    constexpr explicit bezier() = default;

    constexpr explicit bezier(const points_array_t& set) { update(set); }

    [[nodiscard]] point_t constexpr value(coord_t t) const
    {
        if (t < -1 || t > 1) {
            return {0, 0};
        }

        if constexpr (degree_ == degree::cubic) {
            double t2 = t * t;
            double t3 = t2 * t;

            coord_t resx = kx_[3] * t3 + kx_[2] * t2 + kx_[1] * t + kx_[0];
            coord_t resy = ky_[3] * t3 + ky_[2] * t2 + ky_[1] * t + ky_[0];

            return {resx, resy};
        }
    }

    void constexpr update(const points_array_t& set)
    {
        if constexpr (degree_ == degree::cubic) {
            kx_[0] = set[0].x;
            kx_[1] = -set[0].x + 3 * set[1].x;
            kx_[2] = 3 * set[0].x - 6 * set[1].x + 3 * set[2].x;
            kx_[3] = -set[0].x + 3 * set[1].x - 3 * set[2].x + set[3].x;

            ky_[0] = set[0].y;
            ky_[1] = -set[0].y + 3 * set[1].y;
            ky_[2] = 3 * set[0].y - 6 * set[1].y + 3 * set[2].y;
            ky_[3] = -set[0].y + 3 * set[1].y - 3 * set[2].y + set[3].y;
        }
        length_steps_ = 0;
    }

    coord_t
    get_length(const std::uint16_t steps = 100)
    {
        if (steps == 0) {
            return 0;
        }

        if (steps == length_steps_) {
            return length_;
        }

        length_steps_ = steps;
        length_       = 0;

        auto p0 = value(0);

        for (std::uint16_t i{1}; i <= length_steps_; i++) {
            auto    p1  = value(i / static_cast<float>(length_steps_));
            coord_t len = std::sqrt((p1.x - p0.x) * (p1.x - p0.x) + (p1.y - p0.y) * (p1.y - p0.y));
            length_ += len;
            p0 = p1;
        }

        return length_;
    }

    ~bezier() = default;

private:
    std::array<double, static_cast<std::uint8_t>(degree_)> kx_{}, ky_{};

    std::uint16_t length_steps_{};
    coord_t       length_{};
};

using bezier_cubic_f     = bezier<degree::cubic, float>;
using bezier_cubic_d     = bezier<degree::cubic, double>;
using bezier_cubic_int32 = bezier<degree::cubic, std::int32_t>;

}    // namespace loveka::components::math::bezier
