# ECG Detectors

A collection of 7 ECG heartbeat detection algorithms implemented in Python. Developed in conjunction with a new ECG database: http://researchdata.gla.ac.uk/716/. This repository also contains a testing class for the MITDB and the new University of Glasgow database. In addition the module `hrv` provides tools to
analyse heartrate variability.


## Installation

### Linux / Mac

```
pip3 install py-ecg-detectors [--user]
```

### Windows

```
pip install py-ecg-detectors [--user]
```

### From source

```
python3 setup.py install [--user]
```

Use the option `--user` if you don't have system-wise write permission.


## ECG Detector Class Usage

Before the detectors can be used the class must first be initalised with the sampling rate of the ECG recording:

```
from ecgdetectors import Detectors
detectors = Detectors(fs)
```

See `usage_example.py` for an example of how to use the detectors.


## Authors

Luis Howell, luisbhowell@gmail.com

Bernd Porr, bernd.porr@glasgow.ac.uk


## citation / DOI

DOI: 10.5281/zenodo.3353396

https://doi.org/10.5281/zenodo.3353396
