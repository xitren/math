#pragma once

#include <algorithm>
#include <cstdint>

namespace xitren::math {

template <class Type>
static inline Type const&
branchless_select(int const compare, Type const& a, Type const& b)
{
    auto const  ptr_a = reinterpret_cast<std::uintptr_t>(&a);
    auto const  ptr_b = reinterpret_cast<std::uintptr_t>(&b);
    auto        ptr_c = ((compare - 1) & (ptr_b ^ ptr_a)) ^ ptr_a;
    Type const& c     = *(reinterpret_cast<Type const*>(ptr_c));
    return c;
}
}    // namespace xitren::math