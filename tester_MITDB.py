import numpy as np
import pandas as pd
import _tester_utils
import wfdb
import pathlib
import os
from ecgdetectors import Detectors


class MITDB_test:
    """
    Benchmarks detectors with against the MITDB.
    You need to download both the MIT arrhythmia database from: https://alpha.physionet.org/content/mitdb/1.0.0/
    and needs to be placed below this directory: "../mit-bih-arrhythmia-database-1.0.0/"
    """
    def __init__(self):
        current_dir = pathlib.Path(__file__).resolve()
        self.mitdb_dir = str(pathlib.Path(current_dir).parents[1]/'mit-bih-arrhythmia-database-1.0.0')

    def single_classifier_test(self, detector, tolerance=0):
        max_delay_in_samples = 350 / 5
        dat_files = []
        for file in os.listdir(self.mitdb_dir):
            if file.endswith(".dat"):
                dat_files.append(file)
        
        mit_records = [w.replace(".dat", "") for w in dat_files]
        
        results = np.zeros((len(mit_records), 5), dtype=int)

        i = 0
        for record in mit_records:
            progress = int(i/float(len(mit_records))*100.0)
            print("MITDB progress: %i%%" % progress)

            sig, fields = wfdb.rdsamp(self.mitdb_dir+'/'+record)
            unfiltered_ecg = sig[:, 0]  

            ann = wfdb.rdann(str(self.mitdb_dir+'/'+record), 'atr')    
            anno = _tester_utils.sort_MIT_annotations(ann)    

            r_peaks = detector(unfiltered_ecg)

            delay = _tester_utils.calcMedianDelay(r_peaks, unfiltered_ecg, max_delay_in_samples)

            if delay > 1:

                TP, FP, FN = _tester_utils.evaluate_detector(r_peaks, anno, delay, tol=tolerance)
                TN = len(unfiltered_ecg)-(TP+FP+FN)
                
                results[i, 0] = int(record)    
                results[i, 1] = TP
                results[i, 2] = FP
                results[i, 3] = FN
                results[i, 4] = TN

            i = i+1
        
        return results


    def classifer_test_all(self, tolerance=0):

        output_names = ['TP', 'FP', 'FN', 'TN']

        total_records = 0
        for file in os.listdir(self.mitdb_dir):
            if file.endswith(".dat"):
                total_records = total_records + 1

        total_results = np.zeros((total_records, 4*len(_tester_utils.det_names)), dtype=int)

        counter = 0
        for det_name in _tester_utils.det_names:

            print('\n'+det_name)

            result = self.single_classifier_test(_tester_utils.det_from_name(det_name, 360), tolerance=tolerance)
            index_labels = result[:, 0]
            result = result[:, 1:]

            total_results[:, counter:counter+4] = result

            counter = counter+4  

        col_labels = []

        for det_name in _tester_utils.det_names:
                for output_name in output_names:
                    label = det_name+" "+output_name
                    col_labels.append(label)

        total_results_pd = pd.DataFrame(total_results, index_labels, col_labels, dtype=int)            
        total_results_pd.to_csv('results_MITDB'+'.csv', sep=',')

        return total_results_pd
