
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define uint unsigned

extern double _compf(const double t, const uint s, const double * I, const double * a2, const uint len);

int main(int argc, char * argv[])
{
    double * I, * a2, val;
    uint i;
    I  = (double *) malloc(sizeof(double) * 12);
    a2 = (double *) malloc(sizeof(double) * 12);

    for (i = 0; i < 12; ++i)
    {
        I[i]  = 1.f;
        a2[i] = 2.f;
    }

    val = _compf(0.1f, 6, I, a2, 12);

    printf("I[0]: %.1f\na2[0]: %.1f\nValue: %.8f\n", I[0], a2[0], val);

    return 0;
}
