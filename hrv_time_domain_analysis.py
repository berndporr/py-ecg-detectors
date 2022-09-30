#!/usr/bin/python3
# Performs heartrate variation timedomain analysis
#
# It calculates the normalised RMSSD during sitting
# and maths.
#
# This comparison is then run with
# - ground truth (hand corrected R time stamps)
# - Wavelet detector
# - EngZee detector
#
# Via the commandline argument one can choose
# Einthoven II or the ECG from the Chest strap
#
#
# Install https://github.com/berndporr/ECG-GUDB
# via pip.
#

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys
from hrv import HRV
from ecgdetectors import Detectors

from ecg_gudb_database import GUDb

maths_rr_sd = []
maths_error_rr_sd = []
maths_true_sd = []

sitting_rr_sd = []
sitting_error_rr_sd = []
sitting_true_sd = []

total_subjects = 25
subject = []

if len(sys.argv) < 2:
    print("Specify 'e' for Einthoven or 'v' for chest strap ECG.")
    exit(1)

for i in range(total_subjects):
#for i in range(2):
    print(i)
    sitting_class = GUDb(i, 'sitting')
    sitting_class.filter_data()
    maths_class = GUDb(i, 'maths')
    maths_class.filter_data()

    detectors = Detectors(sitting_class.fs)

    if sitting_class.anno_cs_exists and maths_class.anno_cs_exists:
        subject.append(i)

        hrv_class = HRV(sitting_class.fs)

        if "e" in sys.argv[1]:
            ecg_channel_sitting = sitting_class.einthoven_II
            ecg_channel_maths = maths_class.einthoven_II
        elif "v" in sys.argv[1]:
            ecg_channel_sitting = sitting_class.cs_V2_V1
            ecg_channel_maths = maths_class.cs_V2_V1
        else:
            print("Bad argument. Specify 'e' for Einthoven or 'v' for the Chest strap.")
            exit(1)

        r_peaks = detectors.swt_detector(ecg_channel_sitting)
        sitting_rr_sd.append(hrv_class.RMSSD(r_peaks,True))
        r_peaks = detectors.swt_detector(ecg_channel_maths)
        maths_rr_sd.append(hrv_class.RMSSD(r_peaks,True))

        sitting_error_rr = detectors.engzee_detector(ecg_channel_sitting)
        sitting_error_rr_sd.append(hrv_class.RMSSD(sitting_error_rr,True))

        maths_error_rr = detectors.engzee_detector(ecg_channel_maths)
        maths_error_rr_sd.append(hrv_class.RMSSD(maths_error_rr,True))

        maths_true_rr = maths_class.anno_cs
        maths_true_sd.append(hrv_class.RMSSD(maths_true_rr,True))
        
        sitting_true_rr = sitting_class.anno_cs
        sitting_true_sd.append(hrv_class.RMSSD(sitting_true_rr,True))


subject = np.array(subject)
width = 0.4

fig, ax = plt.subplots()
rects1 = ax.bar(subject+(0*width), sitting_true_sd, width)
rects2 = ax.bar(subject+(1*width), maths_true_sd, width)

ax.set_ylabel('SDNN (s)')
ax.set_xlabel('Subject')
ax.set_ylim([0,0.1])
ax.set_title('HRV for sitting and maths test')
ax.set_xticks(subject + width)
ax.set_xticklabels(subject)
ax.legend((rects1[0], rects2[0]), ('sitting', 'maths' ))

plt.figure()

ymax = 0.25

# now let's do stats with no error

avg_sitting_rr_sd = np.average(sitting_rr_sd)
sd_sitting_rr_sd = np.std(sitting_rr_sd)

avg_maths_rr_sd = np.average(maths_rr_sd)
sd_maths_rr_sd = np.std(maths_rr_sd)

plt.bar(['sitting','maths'],[avg_sitting_rr_sd,avg_maths_rr_sd],yerr=[sd_sitting_rr_sd,sd_maths_rr_sd],align='center', alpha=0.5, ecolor='black', capsize=10)
plt.ylim([0,ymax])
plt.title("WAVELET: Sitting vs Maths")
plt.ylabel('nRMSSD')

# and stats with error

avg_sitting_error_rr_sd = np.average(sitting_error_rr_sd)
sd_sitting_error_rr_sd = np.std(sitting_error_rr_sd)

avg_maths_error_rr_sd = np.average(maths_error_rr_sd)
sd_maths_error_rr_sd = np.std(maths_error_rr_sd)

avg_sitting_true_sd = np.average(sitting_true_sd)
sd_sitting_true_sd = np.std(sitting_true_sd)

avg_maths_true_sd = np.average(maths_true_sd)
sd_maths_true_sd = np.std(maths_true_sd)

plt.figure()

plt.bar(['sitting','maths'],[avg_sitting_error_rr_sd,avg_maths_error_rr_sd],yerr=[sd_sitting_error_rr_sd,sd_maths_error_rr_sd],align='center', alpha=0.5, ecolor='black', capsize=10)
plt.ylim([0,ymax])
plt.title("Engzee DETECTOR: Sitting vs Maths")
plt.ylabel('nRMSSD')

plt.figure()

plt.bar(['sitting','maths'],[avg_sitting_true_sd,avg_maths_true_sd],yerr=[sd_sitting_true_sd,sd_maths_true_sd],align='center', alpha=0.5, ecolor='black', capsize=10)
plt.ylim([0,ymax])
plt.title("GROUND TRUTH: Sitting vs Maths")
plt.ylabel('nRMSSD')

t,p = stats.wilcoxon(sitting_true_sd,maths_true_sd)
print("GROUND TRUTH (sitting vs maths): p=",p)

t,p = stats.wilcoxon(sitting_rr_sd,maths_rr_sd)
print("WAVELET (sitting vs maths): p=",p)

t,p = stats.wilcoxon(sitting_error_rr_sd,maths_error_rr_sd)
print("EngZee DETECTOR: (sitting vs maths): p=",p)

t,p = stats.wilcoxon(sitting_true_sd,sitting_rr_sd)
print("Sitting: Wavelet vs ground truth, p=",p)

t,p = stats.wilcoxon(sitting_true_sd,sitting_error_rr_sd)
print("Sitting: EngZee vs ground truth, p=",p)

t,p = stats.wilcoxon(maths_true_sd,maths_rr_sd)
print("Maths: Wavelet vs ground truth, p=",p)

t,p = stats.wilcoxon(maths_true_sd,maths_error_rr_sd)
print("Maths: EngZee vs ground truth, p=",p)

plt.show()
