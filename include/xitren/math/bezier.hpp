#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <utility>

namespace xitren::math {

/**
 * @brief An enum class that represents the degree of a Bezier curve.
 *
 * The degree of a Bezier curve is the number of control points it uses to define the curve.
 * A Bezier curve of degree n has n + 1 control points.
 *
 * The `degree` enum class is an alternative to using an integer to represent the degree of a Bezier curve,
 * as it provides compile-time type safety and prevents invalid degrees from being used.
 *
 * The `degree` enum class has two enumerators: `quadratic` (degree 3) and `cubic` (degree 4).
 */
enum class degree { quadratic = 3, cubic = 4 };

/**
 * @brief A concept that checks if a given type is a valid coordinate type.
 *
 * A coordinate type is a floating-point type (float, double) or an integer type (int, long, etc.)
 *
 * @tparam T The type to be checked.
 * @return `true` if the type is a valid coordinate type, `false` otherwise.
 */
template <typename T>
concept coordinate_typename = std::same_as<float, T> || std::same_as<double, T> || std::same_as<std::int32_t, T>;

/**
 * @brief A class that represents a Bezier curve of a specific degree.
 *
 * The Bezier curve is defined by a set of control points, which define the curve's shape.
 * The degree of the Bezier curve is the number of control points it uses to define the curve.
 *
 * The `bezier` class provides methods for calculating the curve's value at a given position,
 * as well as its length.
 *
 * @tparam degree_ The degree of the Bezier curve. Must be either `degree::quadratic` (degree 3) or `degree::cubic`
 * (degree 4).
 * @tparam T The type of the coordinates of the Bezier curve. Must be a floating-point type (float, double) or an
 * integer type (int, long, etc.).
 */
template <degree degree_, coordinate_typename T>
class bezier {
public:
    static_assert(degree_ == degree::cubic, "Only cubic bezier is supported yet");

    using coord_t = T;

    /**
     * @brief A struct that represents a single control point of the Bezier curve.
     *
     * A control point consists of a `x` and `y` coordinate.
     */
    struct point_t {
        coord_t x{};
        coord_t y{};
    };

    /**
     * @brief An array of control points that define the Bezier curve.
     *
     * The array has a fixed size of `degree_`, which determines the degree of the Bezier curve.
     */
    using points_array_t = std::array<point_t, static_cast<std::uint8_t>(degree_)>;

    /**
     * @brief Default constructor.
     *
     * Creates an empty Bezier curve, with no control points.
     */
    constexpr explicit bezier() = default;

    /**
     * @brief Constructor that initializes the Bezier curve with a set of control points.
     *
     * @param set The set of control points that define the Bezier curve.
     */
    constexpr explicit bezier(points_array_t const& set) { update(set); }

    /**
     * @brief Calculates the value of the Bezier curve at a given position.
     *
     * @param t The position along the curve, where `0` represents the start of the curve and `1` represents the end.
     * @return The value of the curve at the given position.
     */
    [[nodiscard]] constexpr point_t
    value(coord_t t) const
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

    /**
     * @brief Updates the control points of the Bezier curve.
     *
     * @param set The new set of control points that define the Bezier curve.
     */
    constexpr void
    update(points_array_t const& set)
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

    /**
     * @brief Calculates the length of the Bezier curve.
     *
     * @param steps The number of steps to use in the calculation. A higher number of steps will result in a more
     * accurate length calculation, but may be slower.
     * @return The length of the Bezier curve.
     */
    coord_t
    get_length(std::uint16_t const steps = 100)
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

}    // namespace xitren::math
