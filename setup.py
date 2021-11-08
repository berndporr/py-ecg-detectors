#!/usr/bin/python3

from setuptools import setup

with open('README.rst') as f:
    long_description = f.read()

setup(
    name='py-ecg-detectors',
    version='1.1.0',
    description="Seven ECG heartbeat detection algorithms and heartrate variability analysis",
    long_description=long_description,
    author='Luis Howell, Bernd Porr',
    author_email='luisbhowell@gmail.com, bernd.porr@glasgow.ac.uk',
    py_modules=['ecgdetectors','hrv'],
    install_requires=['numpy',
                      'pathlib2;python_version<"3.4"',
                      'scipy',
                      'gatspy',
                      'pywavelets'],
    zip_safe=False,
    url='https://github.com/berndporr/py-ecg-qrs-detectors',
    license='GPL 3.0',
    classifiers=[
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
    ],
)
