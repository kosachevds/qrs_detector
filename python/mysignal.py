class Signal:
    def __init__(self, values, sampling_rate, peaks=None):
        self.values = values
        self.sampling_rate = sampling_rate
        if peaks is None:
            peaks = []
        self.peaks = peaks


def read(filename, begin_sec=0, end_sec=None):
    if end_sec is not None and end_sec <= begin_sec:
        raise ValueError("Invalid end_sec value")
    with open(filename) as in_file:
        sampling_rate = int(in_file.readline().strip())
        begin_line = begin_sec * sampling_rate
        lines = in_file.readlines()[begin_line:]
        if end_sec is not None:
            end_line = end_sec * sampling_rate - begin_line
            lines = lines[:end_line]
    signal = []
    peaks = []
    for index, line in enumerate(lines):
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) == 2:
            is_peak = bool(int(parts[1]))
            if is_peak:
                peaks.append(index)
        signal.append(float(parts[0]))
    return Signal(signal, sampling_rate, peaks)
