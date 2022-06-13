#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import pathlib
from ecgdetectors import Detectors
import sys

current_dir = pathlib.Path(__file__).resolve()

example_dir = current_dir.parent/'example_data'/'ECG.tsv'
unfiltered_ecg_dat = np.loadtxt(example_dir) 
unfiltered_ecg = unfiltered_ecg_dat[:, 0]
fs = 250

detectors = Detectors(fs)

# selected detector by the user (default is the two average one)
seldet = -1

if len(sys.argv) > 1:
    seldet = int(sys.argv[1])
else:
    print("Select another detector by specifying the index as: {} <index>".format(sys.argv[0]))
    print("The following detectors are available:")
    for i in range(len(detectors.get_detector_list())):
        print(i,detectors.get_detector_list()[i][0])
    print("The default detector is the Two Average detector.")

if seldet < 0:
    # default detector
    r_peaks = detectors.two_average_detector(unfiltered_ecg)
else:
    # We use the input argument to select a detector
    r_peaks = detectors.get_detector_list()[seldet][1](unfiltered_ecg)

# If you want to always use the same det then directly call it:
#r_peaks = detectors.two_average_detector(unfiltered_ecg)
#r_peaks = detectors.matched_filter_detector(unfiltered_ecg,"templates/template_250hz.csv")
#r_peaks = detectors.swt_detector(unfiltered_ecg)
#r_peaks = detectors.engzee_detector(unfiltered_ecg)
#r_peaks = detectors.christov_detector(unfiltered_ecg)
#r_peaks = detectors.hamilton_detector(unfiltered_ecg)
#r_peaks = detectors.pan_tompkins_detector(unfiltered_ecg)
#r_peaks = detectors.wqrs_detector(unfiltered_ecg)


plt.figure()
plt.plot(unfiltered_ecg)
plt.plot(r_peaks, unfiltered_ecg[r_peaks], 'ro')
plt.title("Detected R peaks")

plt.show()
