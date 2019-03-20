#include "qrs_detection.h"
#include <stdlib.h>
#include <string.h>

#define WINDOW_SEC (0.160)
#define MIN_RR_SEC (0.200)

static void FilterData(double const* signal, int size, double* output);
static void ComputeDerivative(double const* signal, int size, double* output);
static void ArrayPow2(double* signal, int size);
static void WindowIntegration(double const* signal, int size, double* output, int window_size);
static int Thresholding(const double* integrated, int size, double rate, char* result, int mir_rr_width);
static void Normalize(double* values, int size);

int DetectPeaks(double const* signal, int size, char* result, double rate)
{
    const int WINDOW_SIZE = (int)(WINDOW_SEC * rate);
    const int MIN_RR = (int)(MIN_RR_SEC * rate);
    int count;
    double *buffer;  // filtered signal, integrated signal
    double *derivative;

    buffer = malloc(size * sizeof(double));
    FilterData(signal, size, buffer);
    Normalize(signal, size);

    derivative = malloc(size * sizeof(double));
    ComputeDerivative(buffer, size, derivative);
    Normalize(derivative, size);
    ArrayPow2(derivative, size);

    WindowIntegration(derivative, size, buffer, WINDOW_SIZE);
    free(derivative);
    memset(result, 0, size * sizeof(char));
    count = Thresholding(buffer, size, rate, result, MIN_RR);
    free(buffer);
    return count;
}

/*****************************************************************************/

void FilterData(double const* signal, int size, double* output)
{
    int index;
    double* buffer;

    buffer = malloc(size * sizeof(double));

    /* Low Pass Filter */
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

    /* High Pass Filter */
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
}

void ComputeDerivative(double const* signal, int size, double* output)
{
    int i;

    for (i = 2; i < size - 2; ++i) {
        double value = (-signal[i - 2] - 2 * signal[i - 1] +
                        2 * signal[i + 1] + signal[i + 2]);
        value /= 8.0;
        output[i] = value * value;
    }
    output[0] = output[1] = output[2];
    output[size - 3] = output[size - 2] = output[size - 1];
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
    int index;

    for (index = 0; index < size; ++index) {
        int i;
        double value = 0.0;

        for (i = 0; i < window_size; ++i) {
            if (index - i < 0) {
                break;
            }
            value += signal[index - i];
        }
        output[index] = value / (double)window_size;
    }
}

int Thresholding(const double* integrated, int size, double rate, char* result, int mir_rr_width)
{
    int i, count, previous;
    double peaki, spki, npki, threshold1;

    spki = npki = 0.0;
    count = 0;
    threshold1 = 0.0;
    previous = -1;
    for (i = 1; i < size - 1; ++i) {
        double peaki = integrated[i];
        if (peaki < integrated[i - 1] || peaki <= integrated[i + 1]) {
            continue;
        }

        if (peaki <= threshold1) {
            npki = 0.875 * npki + 0.125 * peaki;
        } else {
            spki = 0.875 * spki + 0.125 * peaki;
        }
        threshold1 = npki + 0.25 * (spki - npki);

        if (peaki > threshold1) {
            if (count == 0) {
                result[i] = 1;
                previous = i;
                ++count;
            } else if (i - previous >= mir_rr_width) {
                result[i] = 1;
                previous = i;
                ++count;
            } else if (integrated[previous] < peaki) {
                result[previous] = 0;
                result[i] = 1;
            }
        }
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
