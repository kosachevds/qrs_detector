from scipy import signal as scisig

_WINDOW_SEC = 0.160
_MIN_RR = 0.2  # compare with 0.33
_ARTICLE_SAMPLING_RATE = 200.0


def detect(signal, rate):
    buffer, samples_delay = _filter_signal(signal, rate)
    buffer = _normalize(buffer)

    buffer = _compute_derivative(buffer)
    buffer = _normalize(buffer)
    buffer = [x * x for x in buffer]

    samples_window = round(_WINDOW_SEC * rate)
    samples_delay += samples_window / 2
    integrated = _window_integration(buffer, samples_window)

    min_rr_samples = round(_MIN_RR * rate)
    indices, _ = _thresholding(integrated, min_rr_samples)
    return [x - samples_delay for x in indices]


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
        # fix: this filters work only for 200 Hz sampling rate
        buffer = _low_pass_filter(signal)
        result = _high_pass_filter(buffer)
        # In the paper delay is 6 samples for LPF and 16 samples for HPF
        # with sampling rate equals 200
        delay = 6 + 16
    else:
        nyq = 0.5 * rate
        lower = 5 / nyq
        upper = 15 / nyq
        b, a = scisig.butter(2, [lower, upper], btype="band")
        result = scisig.filtfilt(b, a, signal)
        delay = 0
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


def _thresholding(integrated, min_rr_samples):
    spki = 0
    npki = 0
    peaks = []
    threshold1 = spki
    th1_list = []
    i = 0
    while i < len(integrated) - 2:
        i += 1
        th1_list.append(threshold1)
        peaki = integrated[i]
        if peaki < integrated[i - 1] or peaki <= integrated[i + 1]:
            continue

        if peaki <= threshold1:
            npki = 0.875 * npki + 0.125 * peaki
        else:
            spki = 0.875 * spki + 0.125 * peaki

        threshold1 = npki + 0.25 * (spki - npki)
        # threshold2 = 0.5 * threshold1

        if peaki > threshold1:
            if not peaks:
                peaks.append(i)
            elif i - peaks[-1] >= min_rr_samples:
                peaks.append(i)
            elif integrated[peaks[-1]] < peaki:
                peaks[-1] = i
    return peaks, th1_list
