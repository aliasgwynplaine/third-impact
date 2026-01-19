#ifndef THIRD_IMPACT_MACROS_H
#define THIRD_IMPACT_MACROS_H

#include <xmmintrin.h>

#define SQUARE(x) ({          \
    __typeof__(x) _x = (x);   \
    _x * _x;                  \
})

#define POW3(x) ({           \
    __typeof__(x) _x = (x);  \
    _x * _x * _x;            \
})

// this can be fast but not so accurate uu
#define FAST_POW(a, b) ({   \
    union {                 \
        double d;           \
        int x[2];           \
    } u = { a };            \
    u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447); \
    u.x[0] = 0;             \
    u.d;                    \
})


#define FAST_RSQRT(x) _mm_cvtss_f32(_mm_rsqrt_ss(_mm_set_ss(x)))

#endif