import os
from matplotlib import pyplot as pp
import detection

_CURRENT_DIR = os.path.dirname(__file__)
_DATA_DIR = os.path.join(_CURRENT_DIR, "..", "data")
_BAD_DATA_FILENAME = "signalwithpeaks.txt"
_SIGNAL_FILENAME = "input"


def main():
    # detect_and_show(_SIGNAL_FILENAME)
    read_signal_with_peaks(_BAD_DATA_FILENAME)


def detect_and_show(filename):
    with open(os.path.join(_DATA_DIR, filename)) as in_file:
        values = in_file.read().split()
        input_signal = [float(x) for x in values]
    result = detection.detect(input_signal, 2000)
    plot_signal_with_peaks(input_signal, result)
    pp.show()


def read_signal_with_peaks(filename):
    with open(os.path.join(_DATA_DIR, filename)) as in_file:
        lines = in_file.readlines()
    signal = []
    peaks = []
    for index, line in enumerate(lines):
        if not line or line == '\n':
            continue
        value, is_peak = line.split()
        signal.append(float(value))
        is_peak = bool(int(is_peak))
        if is_peak:
            peaks.append(index)
    plot_signal_with_peaks(signal, peaks)
    pp.show()


def plot_signal_with_peaks(signal, peaks):
    pp.plot(signal)
    for peak in peaks:
        pp.axvline(peak, color="r")

if __name__ == '__main__':
    main()
