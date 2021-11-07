=============
ECG Detectors
=============

A collection of 7 ECG heartbeat detection algorithms implemented in Python. Developed in conjunction with a new ECG database: http://researchdata.gla.ac.uk/716/. This repository also contains a testing class for the MITDB and the new University of Glasgow database. In addition the module `hrv` provides tools to
analyse heartrate variability.


Installation
============

Linux / Mac::

  pip3 install py-ecg-detectors [--user]

Windows::

  pip install py-ecg-detectors [--user]

From source::

  python3 setup.py install [--user]


Use the option `--user` if you don't have system-wise write permission.


ECG Detector Class Usage
========================

Before the detectors can be used the class must first be initalised with the sampling rate of the ECG recording:

.. code-block:: python

  from ecgdetectors import Detectors
  detectors = Detectors(fs)

See `usage_example.py` for an example of how to use the detectors.

Hamilton
--------

Implementation of P.S. Hamilton, “Open Source ECG Analysis Software Documentation”, E.P.Limited, 2002. Usage::
  
  r_peaks = detectors.hamilton_detector(unfiltered_ecg)

  
Christov
--------

Implementation of Ivaylo I. Christov, “Real time electrocardiogram QRS detection using combined adaptive threshold”, BioMedical Engineering OnLine 2004, vol. 3:28, 2004. Usage::

  r_peaks = detectors.christov_detector(unfiltered_ecg)


Engelse and Zeelenberg
----------------------

Implementation of W. Engelse and C. Zeelenberg, “A single scan algorithm for QRS detection and feature extraction”, IEEE Comp. in Cardiology, vol. 6, pp. 37-42, 1979 with modifications A. Lourenco, H. Silva, P. Leite, R. Lourenco and A. Fred, “Real Time Electrocardiogram Segmentation for Finger Based ECG Biometrics”, BIOSIGNALS 2012, pp. 49-54, 2012. Usage::
  
  r_peaks = detectors.engzee_detector(unfiltered_ecg)



Pan and Tompkins
----------------

Implementation of Jiapu Pan and Willis J. Tompkins. “A Real-Time QRS Detection Algorithm”. In: IEEE Transactions on Biomedical Engineering BME-32.3 (1985), pp. 230–236. Usage::
  
  r_peaks = detectors.pan_tompkins_detector(unfiltered_ecg)


Stationary Wavelet Transform
----------------------------

Implementation based on Vignesh Kalidas and Lakshman Tamil. “Real-time QRS detector using Stationary Wavelet Transform for Automated ECG Analysis”. In: 2017 IEEE 17th International Conference on Bioinformatics and Bioengineering (BIBE). Uses the Pan and Tompkins thresolding method. Usage::
  
  r_peaks = detectors.swt_detector(unfiltered_ecg)


Two Moving Average
------------------

Implementation of Elgendi, Mohamed & Jonkman, Mirjam & De Boer, Friso. (2010). "Frequency Bands Effects on QRS Detection" The 3rd International Conference on Bio-inspired Systems and Signal Processing (BIOSIGNALS2010). 428-431.
Usage::
  
  r_peaks = detectors.two_average_detector(unfiltered_ecg)

  

Matched Filter
--------------

FIR matched filter using template of QRS complex. Uses the Pan and Tompkins thresolding method.
The ECG template is a text file where the samples are in a single column. See
the templates folder on github for examples. Usage::

  r_peaks = detectors.matched_filter_detector(unfiltered_ecg,template_file)

WQRS
--------------

Uses the wqrs detector by Zong, GB Moody, D Jiang. Usage::

  r_peaks = detectors.wqrs_detector(unfiltered_ecg)


Heartrate variability analysis
==============================

The module `hrv` provides a large collection of heartrate
variability measures which are methods of the class `HRV`::

  HR(self, rr_samples)
     Calculate heart-rates from R peak samples.

  NN20(self, rr_samples)
     Calculate NN20, the number of pairs of successive
     NNs that differ by more than 20 ms.

  NN50(self, rr_samples)
     Calculate NN50, the number of pairs of successive
     NNs that differ by more than 50 ms.

  RMSSD(self, rr_samples, normalise=False)
     Calculate RMSSD (root mean square of successive differences).

  SDANN(self, rr_samples, average_period=5.0, normalise=False)
     Calculate SDANN, the standard deviation of the average
     RR intervals calculated over short periods.

  SDNN(self, rr_samples, normalise=False)
     Calculate SDNN, the standard deviation of NN intervals.

  SDSD(self, rr_samples)
     Calculate SDSD (standard deviation of successive differences),
     the standard deviation of the successive differences between adjacent NNs.

  fAnalysis(self, rr_samples)
     Frequency analysis to calc self.lf, self.hf,
     returns the LF/HF-ratio.

  pNN20(self, rr_samples)
     Calculate pNN20, the proportion of NN20 divided by total number of NNs.

  pNN50(self, rr_samples)
     Calculate pNN50, the proportion of NN50 divided by total number of NNs.

For parameters and additional info use the python help function::

  import hrv
  help(hrv)

The example `hrv_time_domain_analysis.py` calculates the heartrate
variability in the timedomain.


Benchmarking
============

`run_all_benchmarks.py` calculates the R peak timestamps
for all detectors, the true/false detections/misses and
saves them in .csv files. Open the script itself or use python's
help function of how to obtain the ECG data such as the MIT db.

`show_stats_plots.py` takes then the .csv files, displays
the results of the different detectors and calculates the stats.

`hrv_time_domain_analysis.py` performs a timedomain analysis
between sitting and a math test using the EngZee detector and
the wavelet detector for comparison.


Realtime / Causal processing
============================
Most ECG R-peak detectors won't detect the actual R-peak so the name
"R-peak detector" is a misnomer. However in practise this won't play
any role as only the temporal differences between R-peaks play a role.
Most detectors work with a threshold which moves the detection forward in time
and use causal filters which delay the detection. Only a
few detectors do actually a maximum detection but even they will be
most likely introducing delays as the ECG will be always filtered by causal
filters. In other words most
detectors cause a delay between the R peak and its detection. That delay
should of course be constant so that the resulting HR and HRV is correct.


Authors
=======

Luis Howell, luisbhowell@gmail.com

Bernd Porr, bernd.porr@glasgow.ac.uk


citation / DOI
==============

DOI: 10.5281/zenodo.3353396

https://doi.org/10.5281/zenodo.3353396
