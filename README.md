# Math header only library for C++20

[![Build passed](https://github.com/xitren/math/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/xitren/math/actions/workflows/cmake-multi-platform.yml)

This project contains a set of math classes for frequent use in embedded projects.

- [Finite Impulse Response (FIR) filter](#fir)
- [Bezier curve](#bezier)
- [Branchless select](#branchless)
- [Kernel-Based Hough Transform](#kht)
- [Matrix Strassen](#strassen)
- [Optimization search](#optimization-search)
- [Proportional Integral Derivative (PID)](#pid-controller)

## Contents

- [Patterns list](#patterns-list)
- [Building and developing](#building-and-developing)
- [Project layout](#project-layout)
- [Contributing](#contributing)
- [Licensing](#licensing)

## Patterns list

### FIR

In signal processing, a finite impulse response (FIR) filter is a filter whose impulse response (or response to any finite length input) is of finite duration, because it settles to zero in finite time. This is in contrast to infinite impulse response (IIR) filters, which may have internal feedback and may continue to respond indefinitely (usually decaying).

Types of FIR Filters

- Low-Pass FIR Filter: It allows frequencies whose cut off frequency is lower than certain value to pass and attenuates those with higher frequencies.

~~~cpp
#include <xitren/math/fir/lowpass.hpp>

fir::lowpass<double, 20, 20, 250>   filter;

filter.value(10.);
~~~

- High-Pass FIR Filter: In contrast, it enables signals having a frequency above a certain cutoff point.

~~~cpp
#include <xitren/math/fir/highpass.hpp>

fir::highpass<double, 20, 20, 250>  filter;

filter.value(10.);
~~~

- Band-Pass FIR Filter: This type of filter allows only some frequencies to go through but prevents others from doing so.

~~~cpp
#include <xitren/math/fir/bandpass.hpp>

fir::bandpass<double, 20, 20, 40, 250> filter;

filter.value(10.);
~~~

- Band-Stop FIR Filter: A filter in this mode lets through all but one narrow range of frequencies.

~~~cpp
#include <xitren/math/fir/bandstop.hpp>

fir::bandstop<double, 20, 20, 40, 250> filter;

filter.value(10.);
~~~

Advantages of FIR Filters

- Stable: Unstable since it does not have a feedback loop.
- Linear Phase: Can be designed to respond linearly in phase.
- Simple Implementation: Can be implemented easily on digital hardware.
- Flexibility: Design can assume various characters of filters.
- No Feedback: Makes the design process for the filter simpler.

Disadvantages of FIR Filters

- High Order: Requires a higher order for sharp frequency responses.
- More Computations: The computations are more intensive due to many coefficients.
- Memory Usage: More memory is used for storing the coefficients.
- Longer Delay: Longer group delay compared to IIR filters.
- Limited Use: Less efficient in real-time processing that requires speed for efficiency reasons.

Applications of FIR Filters

- Audio Signal Processing: These are employed in equalizers and also noise reduction systems.
- Data Transmission: They are vital when designing modems as well as other communication equipment.
- Image Processing: This applies to edge detection, image enhancement, etc.
- Speech Processing: They are used in echo cancellers and voice recognition systems among others.
- Radar Systems: these become important when filtering and detecting signals are concerned.

### Bezier

A Bézier curve is a parametric curve used in computer graphics and related fields. A set of discrete "control points" defines a smooth, continuous curve by means of a formula. Usually the curve is intended to approximate a real-world shape that otherwise has no mathematical representation or whose representation is unknown or too complicated. The Bézier curve is named after French engineer Pierre Bézier (1910–1999), who used it in the 1960s for designing curves for the bodywork of Renault cars.

~~~cpp
#include <xitren/math/bezier.hpp>

bezier_point<int>                p0{-7, 7}, p1{-7, 7}, p2{7, 7}, p3{7, 7};
std::array<bezier_point<int>, 4> base_points{p0, p1, p2, p3};

bezier_quadratic<int, 100> curve{};
curve.update(base_points);
~~~

### Branchless

Branchless programming is a programming technique that eliminates the branches (if, switch, and other conditional statements) from the program.
If you are writing low latency high-performance software, this can be helpful. This can give you a noticeable performance improvement especially if you use conditionals iteratively.
Modern compilers have built-in optimizers which can recognize some patterns and replace them with branch-less counterparts. But these are limited, so it's always better to write optimized code.

~~~cpp
#include <xitren/math/branchless.hpp>

int a = 6;
int b = 7;

auto& sel = branchless_select(a < b, a, b);
std::cout << sel << std::endl;
EXPECT_EQ(sel, a);
~~~

### KHT

Kernel-Based Hough Transform for Detecting Straight Lines in Images.

The KHT is a real-time line detection procedure that extends the conventional voting procedure of the Hough transform. It operates on clusters of approximately collinear pixels. For each cluster, the KHT casts votes using an oriented elliptical-Gaussian kernel that models the uncertainty associated with the best-fitting line for the corresponding cluster. The proposed approach not only significantly improves the performance of the voting scheme, but also produces a much cleaner voting map and makes the transform more robust to the detection of spurious lines.

~~~cpp
#include <xitren/math/kht_opt.hpp>

constexpr std::size_t width  = 16;
constexpr std::size_t height = 16;
vault<width, height>  image{
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0}};

kht<width, height>::convert(image);

image.image();
image.mirror();
~~~

### Strassen

// Later a bit

~~~cpp
~~~

### Optimization search

// Later a bit

~~~cpp
~~~

### PID controller

A proportional–integral–derivative controller (PID controller or three-term controller) is a feedback-based control loop mechanism commonly used to manage machines and processes that require continuous control and automatic adjustment. It is typically used in industrial control systems and various other applications where constant control through modulation is necessary without human intervention.

~~~cpp
#include <xitren/math/pid.hpp>

constexpr double ts{0.1}, ki{0.5};

double setpoint{0.5};
double control_value{};

pid test{ts, 0., 0., ki, 0., 10., -10.};

test.target(setpoint);
~~~

## Building and developing

The recommended way to develop and build this project is to use the Docker image as a dev container.

### Presets

This project makes use of [CMake
presets](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html) to simplify the
process of configuring the project. As a developer, you should use a
`CMakePresets.json` file at the top-level directory of the repository.

### Configure, build and test

You can configure, build and test the `clang_host_release_linux` parts of the project with the following
commands:

~~~shell
cmake --preset=clang_host_release_linux
cmake --build --preset=clang_host_release_linux -t test
~~~

## Project layout

The following ideas are mainly stolen from [P1204R0 – Canonical Project
Structure](https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2018/p1204r0.html).

- Also, all includes, even the "project local" ones use `<>` instead of `""`.
- Subfolders for the source code should comprise somewhat standalone "components".
- There should be a STATIC or INTERFACE library target for each component (this should
  make linking source code dependencies for tests easier).
- Everything test related is in `tests/` and its subdirectories.
- Tests have the `.cpp` extension and are
  named after the class, file, functionality, interface or whatever else they test.
- Hardware tests should be similar to unit tests and check simple functionalities of
  low-level code.
- Golden Tests are used for high level integration/system tests.

The following shows what the directory structure could actually look like.

<details>
  <summary>Directory structure</summary>

  ~~~
  math/
  ├── .github/
  ├── include/
  │   ├── xitren/math/
  │   │   ├── bezier.hpp
  │   │   ├── branchless.hpp
  │   │   ├── kht_opt.hpp
  │   │   ├── matrix_classic.hpp
  │   │   ├── matrix_strassen.hpp
  │   │   ├── optimization.hpp
  │   │   ├── pid.hpp
  │   │   │
  │   │   ├── fir/
  │   │   │   ├── filter.hpp
  │   │   │   ├── lowpass.hpp
  │   │   │   ├── highpass.hpp
  │   │   │   ├── bandstop.hpp
  │   │   │   ├── bandpass.hpp
  │   │   │   └── moving_average.hpp
  │   │   └── ...
  |   └── ...
  │
  ├── tests/
  │   ├── CMakeLists.txt
  │   ├── math_bezier_test.cpp
  │   ├── math_branchless.cpp
  │   ├── math_fir_test.cpp
  │   ├── math_kht_test.cpp
  │   ├── math_matrix_power2_double_test.cpp
  │   ├── math_matrix_power2_int_test.cpp
  │   ├── math_matrix_power2_uint8_test.cpp
  │   ├── math_matrix_test.cpp
  │   ├── math_optimization_test.cpp
  │   ├── math_pid_test.cpp
  │   └── ...
  │
  ├── .clang-format
  ├── .clang-tidy
  ├── .gitignore
  ├── CMakeLists.txt
  ├── CMakePresets.json
  ├── Doxyfile
  ├── LICENSE
  ├── README.md
  └── ...
  ~~~

</details>

## Contributing

The best chance of getting a problem fixed is to submit a patch that fixes it (along with a test case that verifies the fix)!
Feel free to create PR.

## Licensing

See the [LICENSE](LICENSE) document.
