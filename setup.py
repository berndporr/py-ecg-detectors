#!/usr/bin/python3

from setuptools import setup

with open('README.rst') as f:
    long_description = f.read()

setup(
    name='ecgdetectors',
    version='1.0.0',
    description="Seven ECG heartbeat detection algorithms and heartrate variability analysis",
    long_description=long_description,
    author='Luis Howell and Bernd Porr',
    author_email='bernd.porr@glasgow.ac.uk',
    py_modules=['ecgdetectors','hrv'],
    include_package_data=True,
    install_requires=['numpy','pywt','pathlib','scipy','biosppy','gatspy'],
    zip_safe=False,
    url='https://github.com/berndporr/py-ecg-qrs-detectors',
    classifiers=[
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GPL',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Medical',
    ],
)
