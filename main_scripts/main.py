# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 09:50:31 2020

@author: huanyu
"""

from matplotlib import pyplot as plt
import numpy as np
#import scipy
#import time
#import os.path

import myConstant as mc
import myReadfile as mr
from myVariables import (Constant, Global, Lineage)
from my_model_MCMCmultiprocessing import run_model_MCMCmultiprocessing, create_lineage_list_by_pastTag#, add_Bayefactor_2file

MODEL_NAME = mc.MODEL_NAME
LINEAGE_TAG = mc.LINEAGE_TAG
OutputFileDir = mc.OutputFileDir
NUMBER_RAND_NEUTRAL = mc.NUMBER_LINEAGE_MLE

'''
#
# Function randomly chooses "small neutral lineages" and "adpative lineages" from lins
#   
def select_small_lineages(lins, R):
    
    rc = max(10,0.0001*R)#max(5, min(30, 30* 10**(np.log10(R)-5)))
    print('rc=',rc)
    #rc = 30
    lins_choice =[]
    prob_choice = []
    
    for lin in lins:
        if lin.TYPETAG == LINEAGE_TAG['ADP']:
            lins_choice.append(lin)
            log10bf = lin.log10_BayesFactor()
            p = min(1, log10bf)
            prob_choice.append(p)
            #prob_choice.append(0.5)
            
        elif (lin.r0 > 0) and (lin.r0 < rc):
            lins_choice.append(lin)
            prob_choice.append(1-lin.r0/rc)
    
    length = min(NUMBER_RAND_NEUTRAL,len(lins_choice))
    prob_choice = np.asarray(prob_choice)/sum(prob_choice)
    rand_index = np.random.choice(a=len(lins_choice), size=length, replace=False, p=prob_choice)
    lins_ = [ lins_choice[i] for i in list(rand_index)]
    
    return lins_
'''
#
# Function randomly chooses lineages from lins
#   
def select_random_lineages(lins):
    
    lins_choice =[]
    
    for lin in lins:
        if lin.r0 > 0:
            lins_choice.append(lin)
    
    length = min(NUMBER_RAND_NEUTRAL,len(lins_choice))
    rand_index = np.random.choice(a=len(lins_choice), size=length, replace=False)
    lins_ = [ lins_choice[i] for i in list(rand_index)]
    
    return lins_


def run_lineages(lins, start_time, end_time, const, lineage_info):
    #s_bar = []
    if (end_time <= const.T) and (start_time >=0 ):    
        glob = Global(const)    
        
        for current_time in range(start_time, end_time):
            
            if current_time >0 :
                
                # READ LINEAGE FROM THE PAST FILES
                lins = create_lineage_list_by_pastTag(lins, current_time, lineage_info, const)
                
                # UPDATE GLOBAL VARIABLE
                # step1: Choose random lineage for likelihood function
                lins_RAND = select_random_lineages(lins)

                # step2: Maximum likelihood estimate
                glob.UPDATE_GLOBAL(current_time, const, lineage_info, lins_RAND, '2d')

                # run SModel_S for all lineages
                run_dict = {'model_name': MODEL_NAME['SS'], 'lineage_name': lineage_info['lineage_name']}
                run_model_MCMCmultiprocessing(run_dict, lins, glob)
                
    else:
        print(f"the input start_time ={start_time} must >=0 & the end_time ={end_time} must <= total time point {const.T}")
    #print(s_bar)
    #return lins_RAND, glob, s_bar#lins

if __name__ == '__main__':
    
    # ##################################
    # Set your filename and case_name
    # ################################## 
    
    #
    # 1. FileName of Barcode Count data
    # 
    datafilename = 'Data_BarcodeCount_simuMEE_20220213' + '.txt'
    # read file
    lins, totalread, t_cycles = mr.my_readfile(datafilename)
    const = Constant(totalread, t_cycles)

    #
    # 2. Name of This Run Case
    #
    case_name = 'Simulation_20220213'
    #case_name = 'SimuBLT_baseline'
    #case_name = 'SimuBLT_harsh'


    # ##################################
    # Run & output results
    # ##################################

    #
    # 3. Run program
    #
    start_time = 1
    end_time = const.T
    lineage_info =  {'lineage_name': case_name +'_v6'}
    #run_lineages(lins, start_time, end_time, const, lineage_info)

    #
    # 4. Output mean fitness,  Output selection coefficient
    #
    mr.output_global_parameters_BFM(lineage_info,const)
    meanfitness_Bayes_cycle, epsilon_Bayes, t_arr_cycle = mr.read_global_parameters_BFM(lineage_info)

    BFM_result = mr.output_Posterior_parameters_Bayes_v5(lineage_info, datafilename)
    beta = mc.beta
    mr.output_Selection_Coefficient_Bayes_v5(lineage_info, datafilename, BFM_result=BFM_result, beta=beta)
    beta = 4
    mr.output_Selection_Coefficient_Bayes_v5(lineage_info, datafilename, BFM_result=BFM_result, beta=beta)
    beta = 5
    mr.output_Selection_Coefficient_Bayes_v5(lineage_info, datafilename, BFM_result=BFM_result, beta=beta)

    '''
    # read simulation
    bcid_simu_type = {}
    s_simu_adp = {}
    fdirname = '../input/' + 'simuMEE_20220213' + '.npz'
    load = np.load(fdirname)
    s_arr_simu = load['s_arr']
    cell_num_simu = load['cell_num']

    meanfitness_simu_gen = load['meanfitness'] #/ np.log2(mc.D)
    t_simu_gen = [i * np.log2(mc.D) for i in range(len(meanfitness_simu_gen))]
    del load

    for bcid in range(len(s_arr_simu)):
        if s_arr_simu[bcid] > 0:
            bcid_simu_type[bcid] = 'ADP'
        else:
            bcid_simu_type[bcid] = 'NEU'
        s_simu_adp[bcid] = s_arr_simu[bcid]

    bcid_all_result = {}
    for j in range(len(BFM_result)):
        bcid_all_result.update({BFM_result[j]['BCID']: j})
    bcid_adp_fit = []
    for j in range(len(BFM_result)):
        if BFM_result[j]['TYPE'] == 'ADP':
            bcid_adp_fit.append(BFM_result[j]['BCID'])
    # bcid_adp_fit = list(bcid_all_result.keys())
    bcid_adp_simu = [k for k, v in bcid_simu_type.items() if v == 'ADP']

    true_positive_bcid = list(set(bcid_adp_fit).intersection(set(bcid_adp_simu)))
    false_negative_bcid = list(set(bcid_adp_simu) - set(true_positive_bcid))
    false_positive_bcid = list(set(bcid_adp_fit) - set(true_positive_bcid))
    true_negative_bcid = list(set(list(bcid_all_result.keys())) - set(bcid_adp_simu) - set(false_positive_bcid))

    print(len(true_positive_bcid), len(false_negative_bcid), len(false_positive_bcid))

    s_fit_TP = [BFM_result[bcid_all_result[bcid]]['s_mean'] for bcid in true_positive_bcid]
    s_std_fit_TP = [BFM_result[bcid_all_result[bcid]]['s_std'] for bcid in true_positive_bcid]

    s_fit_FP = [BFM_result[bcid_all_result[bcid]]['s_mean'] for bcid in false_positive_bcid]
    s_std_fit_FP = [BFM_result[bcid_all_result[bcid]]['s_std'] for bcid in false_positive_bcid]

    false_negative_bcid_exist = list(set(false_negative_bcid).intersection(set(list(bcid_all_result.keys()))))
    s_fit_FN = [BFM_result[bcid_all_result[bcid]]['s_mean'] for bcid in false_negative_bcid_exist]
    s_std_fit_FN = [BFM_result[bcid_all_result[bcid]]['s_std'] for bcid in false_negative_bcid_exist]

    plt.figure()
    ds = 0.005
    plt.hist(s_fit_TP / np.log2(mc.D), bins=int((max(s_fit_TP) - min(s_fit_TP)) / ds), color='r', alpha=0.2,
             label=f'True Positive (n={len(s_fit_TP)})')
    plt.hist(s_fit_FP / np.log2(mc.D), bins=int((max(s_fit_FP) - min(s_fit_FP)) / ds), color='b', alpha=0.2,
             label=f'False Positive (n={len(s_fit_FP)})')
    plt.hist(s_fit_FN / np.log2(mc.D), bins=int((max(s_fit_FN) - min(s_fit_FN)) / ds), color='m', alpha=0.3,
             label=f'False Negative (n={len(s_fit_FN)})')
    plt.xlabel('s mean (%)')
    plt.legend()
    plt.savefig('../output/DFE_BFM_fit_s_TP_FP_FN_adj' + '.pdf')

    plt.figure()
    ds = 0.005
    plt.hist(s_std_fit_TP / np.log2(mc.D), bins=int((max(s_std_fit_TP) - min(s_std_fit_TP)) / ds), color='r', alpha=0.2,
             label=f'True Positive (n={len(s_std_fit_TP)})')
    plt.hist(s_std_fit_FP / np.log2(mc.D), bins=int((max(s_std_fit_FP) - min(s_std_fit_FP)) / ds), color='b', alpha=0.2,
             label=f'False Positive (n={len(s_std_fit_FP)})')
    plt.hist(s_std_fit_FN / np.log2(mc.D), bins=int((max(s_std_fit_FN) - min(s_std_fit_FN)) / ds), color='m', alpha=0.2,
             label=f'False Positive (n={len(s_std_fit_FN)})')
    plt.xlabel('s std (%)')
    plt.legend()
    plt.savefig('../output/DFE_BFM_fit_s_std_TP_FP_FN' + '.pdf')

    s_fit_All = [BFM_result[j]['s_mean'] for j in range(len(BFM_result))]
    s_std_fit_All = [BFM_result[j]['s_std'] for j in range(len(BFM_result))]
    plt.figure()
    ds = 0.005
    plt.hist(s_fit_All / np.log2(mc.D), bins=int((max(s_fit_All) - min(s_fit_All)) / ds), color='k', alpha=0.2,
             label=f'all(n={len(s_fit_All)})')
    plt.xlabel('s mean (%)')
    plt.legend()
    plt.savefig('../output/DFE_BFM_fit_s_all' + '.pdf')

    plt.figure()
    ds = 0.005
    plt.hist(s_std_fit_All / np.log2(mc.D), bins=int((max(s_std_fit_All) - min(s_std_fit_All)) / ds), color='k',
             alpha=0.2, label=f'all(n={len(s_std_fit_All)})')
    plt.xlabel('s std (%)')
    plt.legend()
    plt.savefig('../output/DFE_BFM_fit_s_std_all' + '.pdf')

    plt.figure()
    plt.plot(s_fit_All / np.log2(mc.D), s_std_fit_All / np.log2(mc.D), 'k.', ms=0.5, alpha=0.02,
             label=f'all (n={len(s_std_fit_All)})')
    plt.plot(s_fit_TP / np.log2(mc.D), s_std_fit_TP / np.log2(mc.D), 'r.', ms=1, alpha=0.3,
             label=f'TP (n={len(s_std_fit_TP)})')
    plt.plot(s_fit_FP / np.log2(mc.D), s_std_fit_FP / np.log2(mc.D), 'b.', ms=1, alpha=0.3,
             label=f'FP (n={len(s_std_fit_FP)})')
    plt.plot(s_fit_FN / np.log2(mc.D), s_std_fit_FN / np.log2(mc.D), 'm.', ms=1, alpha=0.3,
             label=f'FN (n={len(s_std_fit_FN)})')
    plt.plot([0, 0.045], [0, 0.015], ':', color='grey')
    plt.xlim(-0.15, 0.15)
    plt.ylim(0, 0.05)
    plt.xlabel('s mean (%)')
    plt.ylabel('s std (%)')
    plt.legend()
    plt.savefig('../output/BFM_fit_s_mean_2_std' + '.pdf')

    plt.close('all')
    '''

    '''
    make_plot_tmp(false_negative_bcid,'False_Negative', BFM_result, s_simu_adp, cell_num_simu, bcid_all_result, meanfitness_Bayes_cycle)
    make_plot_tmp(false_positive_bcid,'False_Positive', BFM_result, s_simu_adp, cell_num_simu, bcid_all_result, meanfitness_Bayes_cycle)
    make_plot_tmp(ture_positive_bcid,'True_Positive', BFM_result, s_simu_adp, cell_num_simu, bcid_all_result, meanfitness_Bayes_cycle)
    make_plot_tmp(true_negative_bcid[0:100], 'True_Negative', BFM_result, s_simu_adp, cell_num_simu, bcid_all_result, meanfitness_Bayes_cycle)
    '''



