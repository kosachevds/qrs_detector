#include "new_r_detection.h"
#include <stdlib.h>
#include <string.h>

#define WINDOW_SEC (0.150)

typedef struct
{
    int size;
    double rate;
    double const* signal;
    double* filtered;
    double* integrated;
} Detector;

static void FilterData(double const* signal, int size, double* output);
static void SquaredDerivative(double const* signal, int size, double* output);
static void WindowIntegration(double const* signal, int size, double* output, int window_size);
static int Thresholding(Detector const* data, char* result);

int DetectPeaks(double const* signal, int size, char* result, double rate)
{
    int count;
    Detector data = { size, rate, signal, NULL, NULL };
    double* buffer;

    data.filtered = malloc(size * sizeof(double));
    FilterData(signal, size, data.filtered);
    buffer = malloc(size * sizeof(double));
    SquaredDerivative(data.filtered, size, buffer);
    data.integrated = malloc(size * sizeof(double));
    WindowIntegration(buffer, size, data.integrated, (int)(WINDOW_SEC * rate));
    free(buffer);
    memset(result, 0, size * sizeof(char));
    count = Thresholding(&data, result);
    free(data.integrated);
    free(data.filtered);
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

void SquaredDerivative(double const* signal, int size, double* output)
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

int Thresholding(Detector const* data, char* result)
{
    const int MIN_RR = (int)(0.2 * data->rate);
    const int WINDOW = (int)(WINDOW_SEC * data->rate);
    int i, count, last_peak;
    double peaki, spki, npki;

    spki = npki = 0.0;
    peaki = data->integrated[0];
    count = 0;
    last_peak = 0;
    for (i = 1; i < data->size; ++i) {
        double noise, threshold;
        if (data->integrated[i] > peaki) {
            peaki = data->integrated[i];
        }

        noise = data->signal[i] - data->filtered[i];
        npki = (npki * (i - 1) + noise) / (double)i;
        npki = 0.875 * npki + 0.125 * peaki;

        spki = (spki * (i - 1) + data->integrated[i]) / (double)i;
        spki = 0.875 * spki + 0.125 * peaki;

        threshold = 0.5 * (npki + 0.25 * (spki - npki));
        if (data->integrated[i] < threshold || i < WINDOW) {
            continue;
        }
        if (last_peak > 0 && i - last_peak < MIN_RR) {
            continue;
        }
        ++count;
        result[i] = 1;
        last_peak = i;
    }
    return count;
}
