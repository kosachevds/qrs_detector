def detect(signal, rate):
    signal = _low_pass_filter(signal)
    signal = _high_pass_filter(signal)
    signal = _squared_derivative(signal)
    signal = _window_integration(signal, int(0.15 * rate))
    return signal


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
