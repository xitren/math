#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <utility>

namespace xitren::math {

template <typename T>
struct bezier_point {
    T x;
    T y;

    bool
    operator==(bezier_point const& other)
    {
        return (x == other.x) && (y == other.y);
    }
};

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
template <typename T, std::size_t Steps>
class bezier_quadratic : public std::array<bezier_point<T>, Steps> {
public:
    /**
     * @brief An array of control points that define the Bezier curve.
     *
     * The array has a fixed size of `degree_`, which determines the degree of the Bezier curve.
     */
    using points_array_t = std::array<bezier_point<double>, 4>;

    /**
     * @brief Default constructor.
     *
     * Creates an empty Bezier curve, with no control points.
     */
    constexpr bezier_quadratic() = default;

    /**
     * @brief Constructor that initializes the Bezier curve with a set of control points.
     *
     * @param set The set of control points that define the Bezier curve.
     */
    constexpr explicit bezier_quadratic(std::array<bezier_point<T>, 4> const& set) { update(set); }

    /**
     * @brief Updates the control points of the Bezier curve.
     *
     * @param set The new set of control points that define the Bezier curve.
     */
    constexpr void
    update(std::array<bezier_point<T>, 4> const& set)
    {
        k_[0].x = set[0].x;
        k_[1].x = -3 * set[0].x + 3 * set[1].x;
        k_[2].x = 3 * set[0].x - 6 * set[1].x + 3 * set[2].x;
        k_[3].x = -set[0].x + 3 * set[1].x - 3 * set[2].x + set[3].x;

        k_[0].y = set[0].y;
        k_[1].y = -3 * set[0].y + 3 * set[1].y;
        k_[2].y = 3 * set[0].y - 6 * set[1].y + 3 * set[2].y;
        k_[3].y = -set[0].y + 3 * set[1].y - 3 * set[2].y + set[3].y;

        for (int i{}; i < Steps; i++) {
            double t  = static_cast<double>(i + 1) / 100.;
            double t2 = t * t;
            double t3 = t2 * t;

            double resx = k_[3].x * t3 + k_[2].x * t2 + k_[1].x * t + k_[0].x;
            double resy = k_[3].y * t3 + k_[2].y * t2 + k_[1].y * t + k_[0].y;
            auto&  item = this->operator[](i);
            item.x      = static_cast<T>(resx);
            item.y      = static_cast<T>(resy);
        }
    }

    ~bezier_quadratic() = default;

private:
    points_array_t k_{};
};

}    // namespace xitren::math
