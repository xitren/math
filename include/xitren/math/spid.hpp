#pragma once

#include <cstdint>

namespace xitren::math {

template <typename Type>
class s_pid {
public:
    struct params {
        Type kp    = 0;
        Type ki    = 0;
        Type kd    = 0;
        Type limit = 0;
    };

    s_pid(Type kp, Type ki, Type kd, Type limit) : ki_(ki), kp_(kp), kd_(kd), limit_(limit) {}
    s_pid(params par) : s_pid(par.kp, par.ki, par.kd, par.limit) {}
    Type
    calc(Type pid_e)
    {
        Type res = 0;

        // I
        res = i_val + ki_ * pid_e;
        do_limit(res);
        i_val = res;

        // P
        res += kp_ * pid_e;
        do_limit(res);

        // D
        res += kd_ * (pid_e - e_prev);
        do_limit(res);
        e_prev = pid_e;

        return res;
    }
    void
    do_limit(Type& val)
    {
        if (val > limit_)
            val = limit_;
        else if (val < -limit_)
            val = -limit_;
    }

    void
    reset()
    {
        e_prev = 0;
        i_val  = 0;
    }
    ~s_pid() {}

private:
    Type ki_, kp_, kd_;
    Type limit_;

public:
    Type e_prev = 0;
    Type i_val  = 0;
};

}    // namespace loveka::components::math::pid
