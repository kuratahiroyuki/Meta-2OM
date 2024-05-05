#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import openpyxl as px
import pandas as pd
import os

columns_measure= ['Sensitivity', 'Specificity', 'Accuracy', 'MCC', 'AUC', 'Precision', 'Recall', 'F1','AUPRC']

### Input parameter setting ###
outfile_name = "result.xlsx"
deep_method_item = ['LGBM'] # ['RF','SVM','XGB','LGBM']
encode_method_item = ['RCKmer', 'DNC' ,'TNC', 'CKSNAP','ANF','PseEIIP'] #'RCKmer','DNC', 'TNC', 'ENAC', 'binary', 'CKSNAP', 'NCP', 'ANF', 'EIIP' ,'PseEIIP'
################################

for deep_method in deep_method_item :
    val_measure=[]
    test_measure=[]
    for i, encode_method in enumerate(encode_method_item ):

        infile_path = "../data/result/%s/%s" %(deep_method, encode_method)
        infile_name = ["val_measures.csv", "test_measures.csv" ]

        infile1 = infile_path + '/' + infile_name[0] #val
        infile2 = infile_path + '/' + infile_name[1] #test

        val_measure.append(  (pd.read_csv(infile1, index_col=0).iloc[-1].values.tolist())) # means
        test_measure.append( (pd.read_csv(infile2, index_col=0).iloc[-1].values.tolist())) # means

    pd_val_measure  = pd.DataFrame(data=val_measure, index=encode_method_item, columns=columns_measure)
    pd_test_measure = pd.DataFrame(data=test_measure, index=encode_method_item, columns=columns_measure)

    print(pd_val_measure)
    print(pd_test_measure)

    pd_val_test = pd.concat([pd_val_measure, pd_test_measure], axis=0)
   
    if os.path.exists(outfile_name):
        mode_f ='a' #'a' w
    else :
        mode_f ='w'
    with pd.ExcelWriter(outfile_name, engine="openpyxl", mode = mode_f) as writer: 
        pd_val_test.to_excel(writer, sheet_name = deep_method) #index=False, header=False








 


