
#include <math.h>

#ifdef USE_SSE2
#include <emmintrin.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef TILE
#define TILE 2
#endif

#ifndef uint
#define uint unsigned
#endif

#ifdef USE_SSE2

#ifdef _MSC_VER /* visual c++ */
# define ALIGN16_BEG __declspec(align(16))
# define ALIGN16_END 
#else /* gcc or icc */
# define ALIGN16_BEG
# define ALIGN16_END __attribute__((aligned(16)))
#endif

typedef ALIGN16_BEG union {
    double   d[2];
    long     l[2];
    __m128d  v;
} ALIGN16_END vec2;

inline
vec2 pow_pd(const vec2 v, const uint p)
{
    vec2 tot;
    uint i;
    tot.v = _mm_mul_pd(v.v, v.v);
    for (i = 0; i < (p - 2); ++i)
        tot.v = _mm_mul_pd(tot.v, v.v);

    return tot;
}
#endif

double _compf(const double t, const uint s, const double * I, const double * a2, const uint len)
{
    uint i, j, maj;
    double sum = 0., tpi2 = -t * pow(M_PI, 2.);
#ifdef USE_SSE2
    vec2 vi, va, vsum, vpow, vtpi2;
    vsum.v  = _mm_setzero_pd(); 
    vtpi2.v = _mm_load1_pd(&tpi2);
#endif

    // 2. * np.power(pi, (2. * s)) * np.sum(np.multiply(np.multiply(np.power(I, s), a2), np.exp(-I * np.power(np.pi, 2) * t)))

    maj = (len / TILE) * TILE;

    for (i = 0; i < maj; i += TILE)
    {
#ifdef USE_SSE2
        vi.v = _mm_load_pd(&I[i]);
        va.v = _mm_load_pd(&a2[i]);
        vpow = pow_pd(vi, s);
        vsum.v = _mm_add_pd(vsum.v, _mm_mul_pd(vpow.v, __mm_mul_pd(va.v, __vrd2_exp(_mm_mul_pd(vi.v, vtpi2.v)))));
#else
        for (j = 0; j < TILE; ++j)
        {
            sum += pow(I[i+j], s) * a2[i+j] * exp(tpi2 * I[i+j]);
        }
#endif
    } 
    
#ifdef USE_SSE2
    sum += vsum.d[0] + vsum.d[1];
#endif

    for (i = maj; i < len; ++i)
        sum += pow(I[i], s) * a2[i] * exp(tpi2 * I[i]);

    sum *= 2. * pow(M_PI, 2. * s);

    return sum;
}
