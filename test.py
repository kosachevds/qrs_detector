import os
from matplotlib import pyplot as pp
import detection


def main():
    current_dir = os.path.dirname(__file__)
    with open(os.path.join(current_dir, "data", "input")) as in_file:
        values = in_file.read().split()
        input_signal = [float(x) for x in values]
    result = detection.detect(input_signal, 2000)
    pp.plot(input_signal)
    result = [x / 35.0 for x in result]
    pp.plot(result)
    pp.show()

if __name__ == '__main__':
    main()