import numpy as np

_WINDOW_SEC = 0.150
_MIN_RR = 0.2


def detect(signal, rate):
    filtered = _low_pass_filter(signal)  # x1
    filtered = _high_pass_filter(filtered)  # x2
    # x4, x3 - is not squared derivative
    integrated = _squared_derivative(filtered)
    integrated = _window_integration(integrated, int(_WINDOW_SEC * rate))  # x5
    return _thresholding(signal, filtered, integrated, rate)


def _low_pass_filter(signal):
    result = []
    for index, value in enumerate(signal):
        if index >= 1:
            value += 2 * result[index - 1]
        if index >= 2:
            value -= result[index - 2]
        if index >= 6:
            value -= 2 * signal[index - 6]
        if index >= 12:
            value += signal[index - 12]
        result.append(value)
    return result


def _high_pass_filter(signal):
    result = []
    for index, value in enumerate(signal):
        value = -value
        if index >= 1:
            value -= result[index - 1]
        if index >= 16:
            value += 32 * signal[index - 16]
        if index >= 32:
            value += signal[index - 32]
        result.append(value)
    return result


def _squared_derivative(signal):
    result = []
    for index in range(2, len(signal) - 2):
        value = (signal[index + 2] + 2 * signal[index + 1] -
                 signal[index - 2] - 2 * signal[index - 1])
        value /= 8.0
        result.append(value * value)
    return result


def _window_integration(signal, window_size):
    result = []
    for i, _ in enumerate(signal):
        value = 0
        for j in range(window_size):
            if i - j < 0:
                break
            value += signal[i - j]
        result.append(value / float(window_size))
    return result


# TODO: remade
def _thresholding(signal, filtered, integrated, rate):
    peaki = integrated[0]
    spki = 0
    npki = 0
    peaks = [0]
    threshold1 = spki
    for i in range(1, len(integrated)):
        peaki = max(peaki, integrated[i])
        noise = (signal[i] - filtered[i])
        npki = (npki * (i - 1) + noise) / i
        npki = 0.875 * npki + 0.125 * peaki
        spki = (spki * (i - 1) + integrated[i]) / i
        spki = 0.875 * spki + 0.125 * peaki

        threshold1 = npki + 0.25 * (spki - npki)
        threshold2 = 0.5 * threshold1

        if integrated[i] >= threshold2:
            if i - peaks[-1] >= _MIN_RR * rate:
                peaks.append(i)
        #TODO: correct first
    return peaks[1:]
