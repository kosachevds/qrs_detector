#ifndef FILTERING_H
#define FILTERING_H

/* Butterworth filter */
struct Filter;
typedef struct Filter Filter;

typedef enum
{
    FT_HIGH_PASS, FT_LOW_PASS, FT_REJECTION
} FilterType;

#ifdef __cplusplus
extern "C" {
#endif

Filter* InitFilter(double cutoff_hz, double rate, FilterType type);

void FilterData(Filter* filter, double* data, int size);

void CloseFilter(Filter** filter);

void ResetFilter(Filter* filter);

#ifdef __cplusplus
}
#endif

#endif // FILTERING_H
