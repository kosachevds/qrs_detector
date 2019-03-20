#define _USE_MATH_DEFINES
#include "filter.h"
#include <stdlib.h>
#include <math.h>

struct Filter
{
    double b0, b1, b2, a1, a2;
    double x0, x1, x2, y1, y2;
};

static void fill(Filter* filter, double cutoff, int is_hpf);
static void fillRejection(Filter* filter, double cutoff);

/*****************************************************************************/

Filter* InitFilter(double cutoff_hz, double rate, FilterType type)
{
    double cutoff = cutoff_hz / rate;
    Filter* filter = calloc(1, sizeof(Filter));

    if (type == FT_REJECTION) {
        fillRejection(filter, cutoff);
    } else {
        fill(filter, cutoff, (type == FT_HIGH_PASS));
    }
    return filter;
}

void FilterData(Filter* filter, double* data, int size)
{
    int i;

    for (i = 0; i < size; ++i) {
        filter->x2 = filter->x1;
        filter->x1 = filter->x0;
        filter->x0 = data[i];
        data[i] = filter->b0 * filter->x0 + filter->b1 * filter->x1 +
                  filter->b2 * filter->x2 - filter->a1 * filter->y1 -
                  filter->a2 * filter->y2;
        filter->y2 = filter->y1;
        filter->y1 = data[i];
    }
}

void CloseFilter(Filter** filter)
{
    free(*filter);
    *filter = NULL;
}

void ResetFilter(Filter* filter)
{
    filter->x0 = 0;
    filter->x1 = 0;
    filter->x2 = 0;
    filter->y1 = 0;
    filter->y2 = 0;
}

/*****************************************************************************/

void fill(Filter* filter, double cutoff, int is_hpf)
{
    const double B = tan(cutoff * M_PI);
    const double BB = B * B;
    const double S = 1.0 + M_SQRT2 * B + BB;

    if (is_hpf) {
        filter->b0 = 1.0 / S;
        filter->b1 = -2.0 * filter->b0;
    } else {
        filter->b0 = BB / S;
        filter->b1 = 2.0 * filter->b0;
    }
    filter->b2 = filter->b0;
    filter->a1 = 2.0 * (BB - 1.0) / S;
    filter->a2 = (1.0 - M_SQRT2 * B + BB) / S;
}

void fillRejection(Filter* filter, double cutoff)
{
    const double TWO_PI = 2 * M_PI;
    const double MU = 0.005;

    filter->b0 = 1.0 - MU;
    filter->b2 = filter->b0;
    filter->b1 = cos(TWO_PI * cutoff) * (2.0 * MU - 2.0);
    filter->a1 = filter->b1;
    filter->a2 = 1.0 - 2.0 * MU;
}
