import numpy as np
import time
# from matplotlib import pyplot as pp

_WINDOW_SEC = 0.150
_MIN_RR = 0.2


def detect(signal, rate):
    delay = 0
    filtered = _low_pass_filter(signal)
    delay += 6
    filtered = _high_pass_filter(filtered)
    delay += 16
    # pp.figure(0)
    # pp.title("Filtered")
    # pp.plot(filtered)
    squared_der = _squared_derivative(filtered)
    integrated = _window_integration(squared_der, int(_WINDOW_SEC * rate))
    delay += round(_WINDOW_SEC * rate) / 2
    # pp.figure(1)
    # pp.title("Integrated")
    # pp.plot(integrated)
    # pp.figure(2)
    # return _thresholding(signal, filtered, integrated, rate)
    indecies = [x - delay for x in _new_thresholding(integrated, rate)]
    return indecies


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
        first = i - (window_size - 1)
        if first < 0:
            first = 0
        result.append(sum(signal[first:(i + 1)]) / float(window_size))
    return result


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
        # TODO: correct first
    return peaks[1:]


def _new_thresholding(integrated, rate):
    min_interval = int(_MIN_RR * rate)
    peak_indicies = _find_peaks(integrated, limit=0.35, spacing=min_interval)
    # peak_indicies = _find_peaks_(integrated, limit=0.35, spacing=min_interval)
    spki = 0
    npki = 0
    peaks = []
    last_peak = 0
    threshold = 0
    for index in peak_indicies:
        if last_peak > 0 and index - last_peak < min_interval:
            continue
        value = integrated[index]
        if value < threshold:
            npki = 0.875 * npki + 0.125 * value
        else:
            peaks.append(index)
            spki = 0.875 * spki + 0.125 * value
            last_peak = index
        threshold = npki + 0.25 * (spki - npki)
    return peaks


def _find_peaks(data, spacing=1, limit=None):
    """
    Finds peaks in `data` which are of `spacing` width and >=`limit`.
    :param ndarray data: data
    :param float spacing: minimum spacing to the next peak (should be 1 or more)
    :param float limit: peaks should have value greater or equal
    :return array: detected peaks indexes array
    """
    data = np.array(data)
    size = data.size
    x = np.zeros(size + 2 * spacing)
    x[:spacing] = data[0] - 1.e-6
    x[-spacing:] = data[-1] - 1.e-6
    x[spacing:spacing + size] = data
    peak_candidate = np.zeros(size)
    peak_candidate[:] = True
    for s in range(spacing):
        start = spacing - s - 1
        h_b = x[start: start + size]  # before
        start = spacing
        h_c = x[start: start + size]  # central
        start = spacing + s + 1
        h_a = x[start: start + size]  # after
        peak_candidate = np.logical_and(peak_candidate,
                                        np.logical_and(h_c > h_b, h_c > h_a))
    ind = np.argwhere(peak_candidate)
    ind = ind.reshape(ind.size)
    if limit is not None:
        ind = ind[data[ind] > limit]
    return ind


def _find_peaks_(data, spacing, limit):
    size = len(data)
    x = [data[0] - 1.0e-6 for _ in range(spacing)]
    x += [0 for _ in range(size)]
    x += [data[-1] - 1.0e-6 for _ in range(spacing)]
    candidate = [True for _ in range(size)]
    for s in range(spacing):
        start = spacing - s - 1
        h_before = x[start:(start + size)]
        start = spacing
        h_central = x[start:(start + size)]
        start = spacing + s + 1
        h_after = x[start:(start + size)]
        candidate = lists_and(candidate,
                              lists_and(lists_greater(h_central, h_before),
                                        lists_greater(h_central, h_after)))
    return [i for i, x in enumerate(candidate) if x and data[i] > limit]


def lists_and(left, right):
    return [x[0] and x[1] for x in zip(left, right)]


def lists_greater(left, right):
    return [x[0] > x[1] for x in zip(left, right)]