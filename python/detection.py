import math
from scipy import signal as scisig


_WINDOW_SEC = 0.160
_MIN_RR = 0.2  # compare with 0.33
_MAX_RR = 2.0
_ARTICLE_SAMPLING_RATE = 200.0
_LOWER_FILTER_HZ = 5.0
_UPPER_FILTER_HZ = 11.0


def detect(signal, rate):
    buffer, samples_delay = _filter_signal(signal, rate)
    buffer = _normalize(buffer)

    buffer = _compute_derivative(buffer)
    buffer = _normalize(buffer)
    buffer = [x * x for x in buffer]

    samples_window = round(_WINDOW_SEC * rate)
    integrated = _window_integration(buffer, samples_window)
    samples_delay += samples_window // 2

    min_rr_samples = round(_MIN_RR * rate)
    max_rr_samples = round(_MAX_RR * rate)
    indices = _thresholding(integrated, min_rr_samples, max_rr_samples)
    indices = [x - samples_delay for x in indices]
    return _correct_peaks(signal, rate, indices)


def _normalize(values, required_max=1.0):
    max_value = max(values)
    return [item / max_value * required_max for item in values]


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


def _filter_signal(signal, rate):
    result = None
    delay = None
    if rate == _ARTICLE_SAMPLING_RATE:
        buffer = _low_pass_filter(signal)
        result = _high_pass_filter(buffer)
        # In the paper delay is 6 samples for LPF and 16 samples for HPF
        # with sampling rate equals 200
        delay = 6 + 16
    else:
        nyq = 0.5 * rate
        lower = _LOWER_FILTER_HZ / nyq
        upper = _UPPER_FILTER_HZ / nyq
        b, a = scisig.butter(2, upper, btype="low")
        result = scisig.filtfilt(b, a, signal)
        b, a = scisig.butter(2, lower, btype="high")
        result = scisig.filtfilt(b, a, signal)
        delay = int(0.06 * rate)
    return result, delay


def _compute_derivative(signal):
    buffer = []
    max_value = 0.0
    for index in range(2, len(signal) - 2):
        value = (signal[index + 2] + 2 * signal[index + 1] -
                 signal[index - 2] - 2 * signal[index - 1])
        value /= 8.0
        if value > max_value:
            max_value = value
        buffer.append(value)
    return buffer


def _window_integration(signal, window_size):
    result = []
    value = 0
    for i, x in enumerate(signal):
        first = i - (window_size - 1)
        value += x / window_size
        if first > 0:
            value -= signal[first - 1] / window_size
        result.append(value)
    return result


def _thresholding(integrated, min_rr_width, max_rr_width):
    spki = 0
    npki = 0
    peaks = []
    threshold1 = spki
    threshold2 = spki
    searchback = False
    searchback_end = 0
    previous = 0
    i = 2
    while i < len(integrated) - 2:
        if i - previous > max_rr_width and i - searchback_end > max_rr_width:
            searchback = True
            searchback_end = i
            i = previous + 2
            continue
        if searchback and i == searchback_end:
            searchback = False
            continue
        peaki = integrated[i]
        if peaki < integrated[i - 2] or peaki <= integrated[i + 2]:
            i += 1
            continue

        is_qrs = False
        if searchback:
            if peaki > threshold2:
                spki = 0.750 * spki + 0.250 * peaki
                is_qrs = True
        elif peaki > threshold1:
            spki = 0.875 * spki + 0.125 * peaki
            is_qrs = True

        if is_qrs:
            if previous == 0 or i - previous >= min_rr_width:
                peaks.append(i)
            elif integrated[previous] < peaki:
                peaks[-1] = i
            previous = i
        else:
            npki = 0.875 * npki + 0.125 * peaki

        threshold1 = npki + 0.25 * (spki - npki)
        threshold2 = 0.5 * threshold1
        i += 1
    return peaks


def _correct_peaks(signal, rate, peaks):
    left_add = int(0.075 * rate)
    right_add = int(0.075 * rate)
    i = 0
    # TODO: debug
    while i < len(peaks):
        old_index = peaks[i]
        begin = max(old_index - left_add, 1)
        end = min(old_index + right_add, len(signal) - 1)
        baseline = (signal[begin] + signal[end]) / 2
        max_value = math.fabs(signal[old_index] - baseline)
        new_index = old_index
        for j in range(begin, end):
            value = math.fabs(signal[j] - baseline)
            if value > max_value:
                max_value = value
                new_index = j
        if new_index != old_index:
            peaks[i] = new_index
        i += 1
    return peaks
