#!/usr/bin/python3

from setuptools import setup

with open('README.rst') as f:
    long_description = f.read()

setup(
    name='py-ecg-detectors',
    version='0.9.1',
    description="Seven ECG heartbeat detection algorithms and timedomain heartrate variability analysis",
    long_description=long_description,
    author='Luis Howell, Bernd Porr',
    author_email='luisbhowell@gmail.com, bernd.porr@glasgow.ac.uk',
    py_modules=['ecgdetectors','hrv'],
    include_package_data=True,
    install_requires=['numpy','pywt','pathlib','scipy','biosppy','gatspy'],
    zip_safe=False,
    url='https://github.com/berndporr/py-ecg-qrs-detectors',
    license='GPL 3.0',
    classifiers=[
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
    ],
)
