_WINDOW_SEC = 0.160
_MIN_RR = 0.2


def detect(signal, rate):
    filtered = _low_pass_filter(signal)
    filtered = _high_pass_filter(filtered)
    squared_derivative = _squared_derivative(filtered)
    samples_window = round(_WINDOW_SEC * rate)
    integrated = _window_integration(squared_derivative, samples_window)

    # In the paper delay is 6 samples for LPF and 16 samples for HPF
    # with sampling rate equals 200
    delay_sec = (6 + 16) / 2000.0
    # delay_sec += _WINDOW_SEC / 2.0
    offset = round(delay_sec * rate)

    min_rr_samples = round(_MIN_RR * rate)
    # indices = _new_thresholding(integrated, min_rr_samples)
    indices = _thresholding(integrated, min_rr_samples)
    _debug_plotting(signal, integrated, indices, offset)
    return [x - offset for x in indices]


def _debug_plotting(signal, integrated, indices, offset):
    from matplotlib import pyplot as pp

    pp.plot(signal)

    signal_max = max(signal)
    integrated_max = max(integrated)
    normalized_integrated = [item / integrated_max * signal_max for item in integrated]
    pp.plot(normalized_integrated)

    for peak in indices:
        pp.axvline(peak, color="r")

    indices_with_offset = [x - offset for x in indices]
    for peak in indices_with_offset:
        pp.axvline(peak, color="g")
    pp.show()


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
    peaks = [0]
    threshold1 = spki
    for i in range(1, len(integrated) - 1):
        peaki = integrated[i]
        if peaki < integrated[i - 1] or peaki < integrated[i + 1]:
            continue

        if peaki <= threshold1:
            npki = 0.875 * npki + 0.125 * peaki
        else:
            spki = 0.875 * spki + 0.125 * peaki

        threshold1 = npki + 0.25 * (spki - npki)
        # threshold2 = 0.5 * threshold1

        if peaki > threshold1 and i - peaks[-1] >= min_rr_samples:
            peaks.append(i)
        # TODO: correct first
    return peaks[1:]
