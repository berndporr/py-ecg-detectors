"""
This script benchmarks all detectors with both the MIT database and the
GUdatabase. For the GU database it runs it for both the chest strap and
Einthoven II with loose cables. You need to download both the MIT arrhythmia
database from: https://alpha.physionet.org/content/mitdb/1.0.0/ and the GU
database from: http://researchdata.gla.ac.uk/716/ Both need to be placed below
this directory: "../mit-bih-arrhythmia-database-1.0.0/" for the MITDB and and
"../dataset_716" for the GU database.
"""
import numpy as np
from multiprocessing import Process
from tester_MITDB import MITDB_test
from tester_GUDB import GUDB_test
from ecgdetectors import Detectors

# benchmark the detectors with the MIT DB
do_test_MIT = True

# benchmark the detectors with the GU DB
do_test_GU = True


def run_GUDB_tests(leads):
    # GUDB database testing
    gu_test = GUDB_test()
    gu_test.classifer_test_all(tolerance=0, config=leads)


def run_MIT_tests():
    # MIT-BIH database testing
    mit_test = MITDB_test()
    mit_detectors = Detectors(360)

    # test single detector
    matched_filter_mit = mit_test.single_classifier_test(
        mit_detectors.matched_filter_detector, tolerance=0)
    np.savetxt('matched_filter_mit.csv', matched_filter_mit, fmt='%i',
               delimiter=',')

    # test all detectors on MITDB, save results to csv, will take some time
    mit_test.classifer_test_all()


if do_test_MIT:
    pmit = Process(target=run_MIT_tests)
    pmit.start()

if do_test_GU:
    pgustrap = Process(target=run_GUDB_tests, args=('chest_strap',))
    pgustrap.start()
    pgucables = Process(target=run_GUDB_tests, args=('loose_cables',))
    pgucables.start()

if do_test_MIT:
    pmit.join()

if do_test_GU:
    pgustrap.join()
    pgucables.join()
