import os
from matplotlib import pyplot as pp
import detection
import mysignal as sig

_CURRENT_DIR = os.path.dirname(__file__)
_DATA_DIR = os.path.join(_CURRENT_DIR, "..", "data")
_INPUT_FILENAME = os.path.join(_DATA_DIR, "input_1209.txt")


def main():
    show_peaks(_INPUT_FILENAME)


def show_peaks(filename, show_file_peaks=False, end_sec=None):
    signal = sig.read(filename, end_sec=end_sec)
    peaks = detection.detect(signal.values, signal.sampling_rate)
    plot_signal_with_peaks(signal.values, peaks, "r")
    if signal.peaks and show_file_peaks:
        plot_vlines(signal.peaks, "g")
    pp.show()


def plot_signal_with_peaks(signal, peaks, peaks_color=None):
    pp.plot(signal)
    if peaks_color is None:
        plot_vlines(peaks)
    else:
        plot_vlines(peaks, peaks_color)


def plot_vlines(x_list, color="r"):
    for item in x_list:
        pp.axvline(item, color=color)


if __name__ == '__main__':
    main()
