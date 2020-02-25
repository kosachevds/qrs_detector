#include "qrs_detection.h"
#include "filter.h"
#include <stdlib.h>
#include <string.h>

#define WINDOW_SEC (0.160)
#define MIN_RR_SEC (0.200)
#define MAX_RR_SEC (2.0)
#define USE_ARTICLE_FILTER (1)
#define LOWER_HZ (5.0)
#define UPPER_HZ (11.0)
#define SUM_FILTER_DELAY_SEC (0.06)

static int FilterSignal(double const* signal, int size, double rate, double* output);
static int ApplyArticleFilter(double const* signal, int size, double* output);
static void ComputeDerivative(double const* signal, int size, double* output);
static void ArrayPow2(double* signal, int size);
static void WindowIntegration(double const* signal, int size, double* output, int window_size);
static int Thresholding(const double* integrated, int size, int min_rr_width, int max_rr_width, char* result);
static void Normalize(double* values, int size);
static void SubtractDelay(char* qrs_detection_result, int size, int samples_delay);

/*****************************************************************************/

int DetectQrsPeaks(double const* signal, int size, char* result, double rate)
{
    const int WINDOW_SIZE = (int)(WINDOW_SEC * rate);
    const int MIN_RR = (int)(MIN_RR_SEC * rate);
    const int MAX_RR = (int)(MAX_RR_SEC * rate);
    int count, delay;
    double *buffer;  // filtered signal, integrated signal
    double *derivative;

    buffer = malloc(size * sizeof(double));
    delay = FilterSignal(signal, size, rate, buffer);
    Normalize(buffer, size);

    derivative = malloc(size * sizeof(double));
    ComputeDerivative(buffer, size, derivative);
    Normalize(derivative, size);
    ArrayPow2(derivative, size);

    WindowIntegration(derivative, size, buffer, WINDOW_SIZE);
    delay += WINDOW_SIZE / 2;
    free(derivative);
    count = Thresholding(buffer, size, MIN_RR, MAX_RR, result);
    free(buffer);
    SubtractDelay(result, size, delay);
    return count;
}

/*****************************************************************************/

int ApplyArticleFilter(double const* signal, int size, double* output)
{
    int index;
    double* buffer;

    buffer = malloc(size * sizeof(double));

    /* Low Pass Filter ~ 11 Hz */
    for (index = 0; index < size; ++index) {
        double value = signal[index];
        if (index >= 1) {
            value += 2 * buffer[index - 1];
        }
        if (index >= 2) {
            value -= buffer[index - 2];
        }
        if (index >= 6) {
            value -= 2 * signal[index - 6];
        }
        if (index >= 12) {
            value += signal[index - 12];
        }
        buffer[index] = value;
    }

    /* High Pass Filter ~ 5 Hz */
    for (index = 0; index < size; ++index) {
        double value = -buffer[index];
        if (index >= 1) {
            value -= output[index - 1];
        }
        if (index >= 16) {
            value += 32 * buffer[index - 16];
        }
        if (index >= 32) {
            value += buffer[index - 32];
        }
        output[index] = value;
    }

    free(buffer);
    return 6 + 16;
}

int FilterSignal(double const* signal, int size, double rate, double* output)
{
    const double ARTICLE_SAMPLING_RATE = 200.0;
    Filter* filter;

#ifdef USE_ARTICLE_FILTER
    if (rate == ARTICLE_SAMPLING_RATE) {
        return ApplyArticleFilter(signal, size, output);
    }
#endif // USE_ARTICLE_FILTER

    memcpy(output, signal, size * sizeof(double));

    filter = InitFilter(UPPER_HZ, rate, FT_LOW_PASS);
    FilterData(filter, output, size);
    CloseFilter(&filter);

    filter = InitFilter(LOWER_HZ, rate, FT_HIGH_PASS);
    FilterData(filter, output, size);
    CloseFilter(&filter);
    return (int)round(SUM_FILTER_DELAY_SEC * rate);
}

void ComputeDerivative(double const* signal, int size, double* output)
{
    int i;

    for (i = 2; i < size - 2; ++i) {
        double value = (-signal[i - 2] - 2 * signal[i - 1] +
                        2 * signal[i + 1] + signal[i + 2]);
        output[i] = value / 8.0;
    }
    output[0] = output[1] = output[2];
    output[size - 1] = output[size - 2] = output[size - 3];
}

void ArrayPow2(double* signal, int size)
{
    int i;

    for (i = 0; i < size; ++i) {
        signal[i] *= signal[i];
    }
}

void WindowIntegration(double const* signal, int size, double* output, int window_size)
{
    int i;
    double value = 0.0;
    //double inv_window_size = 1.0 / window_size;

    for (i = 0; i < size; ++i) {
        int first = i - (window_size - 1);
        //value += signal[index] * inv_window_size;
        value += signal[i] / window_size;
        if (first > 0) {
            //value -= signal[first - 1] * inv_window_size;
            value -= signal[first - 1] / window_size;
        }
        output[i] = value;
    }
}

int Thresholding(const double* integrated, int size, int min_rr_width, int max_rr_width, char* result)
{
    int i, count, previous, searchback, searchback_end;
    double spki, npki, threshold1, threshold2;

    spki = npki = 0.0;
    count = 0;
    threshold1 = 0.0;
    threshold2 = 0.0;
    previous = 0;
    searchback_end = 0;
    searchback = 0;
    memset(result, 0, size * sizeof(char));
    for (i = 2; i < size - 2; ++i) {
        int is_qrs;
        double peaki;
        if (i - previous > max_rr_width && i - searchback_end > max_rr_width) {
            searchback = 1;
            searchback_end = i;
            i = previous + 1;
            continue;
        }
        if (searchback && i == searchback_end) {
            searchback = 0;
            continue;
        }
        peaki = integrated[i];
        if (peaki < integrated[i - 2] || peaki <= integrated[i + 2]) {
            continue;
        }
        is_qrs = 0;
        if (searchback) {
            if (peaki > threshold2) {
                spki = 0.750 * spki + 0.250 * peaki;
                is_qrs = 1;
            }
        } else if (peaki > threshold1) {
            spki = 0.875 * spki + 0.125 * peaki;
            is_qrs = 1;
        }
        if (is_qrs) {
            if (count == 0 || i - previous >= min_rr_width) {
                result[i] = MARK_QRS;
                ++count;
            } else if (integrated[previous] < peaki) {
                result[previous] = MARK_NO_QRS;
                result[i] = MARK_QRS;
            }
            previous = i;
        } else {
            npki = 0.875 * npki + 0.125 * peaki;
        }
        threshold1 = npki + 0.25 * (spki - npki);
        threshold2 = 0.5 * threshold1;
        ++i;
    }
    return count;
}

void Normalize(double* values, int size)
{
    int i;
    double max_value = values[0];

    for (i = 1; i < size; ++i) {
        if (values[i] > max_value) {
            max_value = values[i];
        }
    }

    for (i = 0; i < size; ++i) {
        values[i] /= max_value;
    }
}

void SubtractDelay(char* qrs_detection_result, int size, int samples_delay)
{
    int i;

    for (i = samples_delay; i < size; ++i) {
        if (qrs_detection_result[i] == MARK_NO_QRS) {
            continue;
        }
        qrs_detection_result[i] = MARK_NO_QRS;
        qrs_detection_result[i - samples_delay] = MARK_QRS;
    }
}
