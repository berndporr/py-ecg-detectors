#!/usr/bin/python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import _tester_utils

def str_join(a, b, c):

    a = np.array(a)
    b = np.array(b)
    c = np.array(c)

    str_array = np.core.defchararray.add(np.core.defchararray.add(a, b), c)

    return str_array


def get_sensivities(results_df, detector_name, experiment=None):

    if experiment!=None:
        tp_col_names = str_join(detector_name+' ', [experiment], ' TP')
        fp_col_names = str_join(detector_name+' ', [experiment], ' FP')
        fn_col_names = str_join(detector_name+' ', [experiment], ' FN')
        tn_col_names = str_join(detector_name+' ', [experiment], ' TN')
        total_tp = (results_df.loc[:, tp_col_names].values)[:,0]
        total_fn = (results_df.loc[:, fn_col_names].values)[:,0]
    else:
        tp_col_names = detector_name+' '+'TP'
        fp_col_names = detector_name+' '+'FP'
        fn_col_names = detector_name+' '+'FN'
        tn_col_names = detector_name+' '+'TN'      
        total_tp = results_df.loc[:, tp_col_names].values
        total_fn = results_df.loc[:, fn_col_names].values

    se = []

    for tp,fn in zip(total_tp,total_fn):
            if (tp + fn) > 0:
                s = tp/(tp+fn)*100.0
                se.append(s)

    return np.array(se)


def get_result(results_df, det_names, experiment=None):
    
    m = []
    s = []
    for det in det_names:
        m.append(np.mean(get_sensivities(results_df, det, experiment)))
        s.append(np.std(get_sensivities(results_df, det, experiment)))

    return np.array(m),np.array(s)



def compare_det_test(results_df, detector_name1, detector_name2, experiment=None):
    se1 = get_sensivities(results_df, detector_name1, experiment)
    if len(se1) < 2:
        return 0
    se2 = get_sensivities(results_df, detector_name2, experiment)
    if len(se2) < 2:
        return 0
    l = min(len(se1),len(se2))
    if (l < 20):
        return None
    #print("1:",se1[:l])
    #print("2:",se2[:l])
    try:
        t,p = stats.wilcoxon(se1[:l],se2[:l])
        return p
    except:
        return None

def compare_cables(results_df1,results_df2,det1,experiment):
    se1 = get_sensivities(results_df1, det1, experiment)
    if len(se1) < 2:
        return 0
    se2 = get_sensivities(results_df2, det1, experiment)
    if len(se2) < 2:
        return 0
    l = min(len(se1),len(se2))
    if (l < 20):
        return None
    try:
        t,p = stats.wilcoxon(se1[:l],se2[:l])
        return p
    except:
        return None


def compare_experiments(results_df1,det1):
    se1 = get_sensivities(results_df1, det1, 'sitting')
    if len(se1) < 2:
        return 0
    se2 = get_sensivities(results_df1, det1, 'jogging')
    if len(se2) < 2:
        return 0
    l = min(len(se1),len(se2))
    if (l < 20):
        return None
    try:
        t,p = stats.wilcoxon(se1[:l],se2[:l])
        return p
    except:
        return None


def single_plot(data, std, y_label, title = None):
    fig, ax = plt.subplots()
    _tester_utils.det_fullnames = ['Elgendi et al', 'Matched Filter', 'Kalidas and Tamil', 'Engzee Mod', 'Christov', 'Hamilton', 'Pan and Tompkins']
    x_pos = np.arange(len(_tester_utils.det_fullnames))

    fig.set_size_inches(10, 7)
    rects1 = ax.bar(x_pos, data, yerr=std, width = 0.65, align='center', alpha=0.5, ecolor='black', capsize=10)
    ax.set_ylim([0,100])
    ax.set_ylabel(y_label)
    ax.set_xlabel('Detector')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(_tester_utils.det_fullnames)

    if title!=None:
        ax.set_title(title)

    plt.tight_layout()

    return rects1


def double_plot(data1, std1, data2, std2, y_label, legend1, legend2, title=None):
    fig, ax = plt.subplots()
    _tester_utils.det_fullnames = ['Elgendi et al', 'Matched Filter', 'Kalidas and Tamil', 'Engzee Mod', 'Christov', 'Hamilton', 'Pan and Tompkins']
    x_pos = np.arange(len(_tester_utils.det_fullnames))

    fig.set_size_inches(10, 7)
    width = 0.4
    rects1 = ax.bar(x_pos, data1, width, yerr=std1, alpha=0.5, ecolor='black', capsize=10)
    rects2 = ax.bar(x_pos+width, data2, width, yerr=std2, alpha=0.5, ecolor='black', capsize=10)
    ax.set_ylim([0,100])
    ax.set_ylabel(y_label)
    ax.set_xlabel('Detector')
    ax.set_xticks(x_pos + width / 2)
    ax.set_xticklabels(_tester_utils.det_fullnames)
    ax.legend((rects1[0], rects2[0]), (legend1, legend2))

    if title!=None:
        ax.set_title(title)

    plt.tight_layout()

    return rects1, rects2


def print_stat(p):
    if p == None:
        print('--- & ',end='')
        return
    s = ""
    if p < 0.05:
        s = "*"
    print('{:03.2f}{} & '.format(p,s),end='')

def print_detector_comparision_table(title,datafile,experiment=None):
    print(title)
    print("      & ",end='')
    for detlong1 in _tester_utils.det_fullnames:
        print(detlong1," & ",end='')
    print("\\\\")
    for det1,detlong1 in zip(_tester_utils.det_names,_tester_utils.det_fullnames):
        print(detlong1," & ",end='')
        for det2,detlong2 in zip(_tester_utils.det_names,_tester_utils.det_fullnames):
            p = compare_det_test(datafile, det1, det2, experiment)
            print_stat(p)
        print("\\\\")
    print()

# GUDB
gudb_cs_results = pd.read_csv('results_GUDB_chest_strap.csv', dtype=int, index_col=0)
gudb_cable_results = pd.read_csv('results_GUDB_loose_cables.csv', dtype=int, index_col=0)

# MITDB
mitdb_results = pd.read_csv('results_MITDB.csv', dtype=int, index_col=0)

experiment_names = ['sitting','maths','walking','hand_bike','jogging']
output_names = ['TP', 'FP', 'FN', 'TN']

print_detector_comparision_table("MIT database",mitdb_results)

print_detector_comparision_table("EINTHOVEN SITTING",gudb_cable_results,'sitting')
print_detector_comparision_table("CHEST STRAP SITTING",gudb_cs_results,'sitting')

print_detector_comparision_table("EINTHOVEN JOGGING",gudb_cable_results,'jogging')
print_detector_comparision_table("CHEST STRAP JOGGING",gudb_cs_results,'jogging')

print("\n\n\n")

print("CHEST STRAP VS EINTH")
for exp1 in ['sitting','jogging']:
    print(exp1,":")
    for det1 in _tester_utils.det_names:
        p = compare_cables(gudb_cable_results,gudb_cs_results,det1,exp1)
        print_stat(p)
    print("\\\\")

print("\n\n\n")

print("LOOSE CABLES: SITTING VS JOGGING")
for det1 in _tester_utils.det_names:
    p = compare_experiments(gudb_cable_results,det1)
    print_stat(p)
print("\\\\")

print("\n\n\n")

print("CHEST STRAP: SITTING VS JOGGING")
for det1 in _tester_utils.det_names:
    p = compare_experiments(gudb_cs_results,det1)
    print_stat(p)
print("\\\\")

print("\n\n\n")

# calculating results
mitdb_avg,mitdb_std = get_result(mitdb_results, _tester_utils.det_names)
print("mitdb:",mitdb_avg)
gudb_cs_sitting_avg,gudb_cs_sitting_std = get_result(gudb_cs_results,
                                                     _tester_utils.det_names, 'sitting')
print("chest strap sitting:",gudb_cs_sitting_avg)
gudb_cable_sitting_avg,gudb_cable_sitting_std = get_result(gudb_cable_results,
                                                           _tester_utils.det_names, 'sitting')
print("lose cables sitting:",gudb_cable_sitting_avg)
gudb_cs_jogging_avg,gudb_cs_jogging_std = get_result(gudb_cs_results,
                                                     _tester_utils.det_names, 'jogging')
gudb_cable_jogging_avg,gudb_cable_jogging_std = get_result(gudb_cable_results,
                                                           _tester_utils.det_names, 'jogging')

# plotting
single_plot(mitdb_avg, mitdb_std, 'Sensitivity (%)', 'MITDB')

double_plot(gudb_cs_sitting_avg, gudb_cs_sitting_std,
            gudb_cable_sitting_avg, gudb_cable_sitting_std,
            'Sensitivity (%)', 'Chest Strap', 'Loose Cables', 'GUDB: cable, sitting')

double_plot(gudb_cs_jogging_avg, gudb_cs_jogging_std,
            gudb_cable_jogging_avg, gudb_cable_jogging_std,
            'Sensitivity (%)', 'Chest Strap', 'Loose Cables', 'GUDB: cable, jogging')



plt.show()
