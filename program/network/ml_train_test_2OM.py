#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 2 09:22:16 2023
@author: Kurata Laboratory
"""

import os
import time
import argparse
import pandas as pd
import numpy as np
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn import svm
import xgboost as xgb
import lightgbm as lgb
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from catboost import CatBoostClassifier
from gensim.models import word2vec
#from libml.encodingRNA_1 import PseEIIP, CTD, NPS, NPPS, NCP_ND, BPB, PseKNC
from encodingRNA_1 import PseEIIP, CTD, NPS, NPPS, NCP_ND, BPB, PseKNC
from encodingRNA_2 import Kmer, RCKmer, DNC, TNC, ENAC, binary, CKSNAP, NCP, ANF, EIIP, PseEIIP, PSTNPss  

def emb_seq_w2v(seq_mat, w2v_model, num): ###
    num_sample = len(seq_mat)    
    seq_emb=[]  
    for j in range(num_sample):
      if j%1000==0:
        print(f'j: {j}')
      seq=seq_mat[j]
      seq_emb.append( np.array([np.array(w2v_model.wv[seq[i:i+num]]) for i in range(len(seq) - num + 1)]) )
    seq_emb = np.array(seq_emb)        
    seq_emb = seq_emb.reshape(num_sample,len(seq) - num + 1, -1)  #74446 x 41 x 4
    seq_emb = seq_emb.reshape(num_sample,1,-1).squeeze() #74446 x 41 x 4
    return seq_emb


def pad_input_csv(filename, seqwin, index_col = None):
    df1 = pd.read_csv(filename, delimiter=',',index_col = index_col)
    seq = df1.loc[:,'seq'].tolist()
    #data triming and padding
    for i in range(len(seq)):
       if len(seq[i]) > seqwin:
         seq[i]=seq[i][0:seqwin]
       seq[i] = seq[i].ljust(seqwin, '-')
    df1['seq'] = seq

    return df1


def pad_input_csv_rnafm(filename, seqwin, index_col = None):
    df1 = pd.read_csv(filename, delimiter=',',index_col = index_col)
    seqs = df1.loc[:,'seq'].tolist()
    #data triming and padding
    sequence = []
    rna_fm_token = []
    for seq in seqs:
        if len(seq) > seqwin:
            seq = seq[0:seqwin]
        sequence.append( seq.ljust(seqwin, '-') )        
        rna_fm_token.append([rna_dict[res] for res in seq])
      
    df1['seq'] = sequence
    df1['token'] = rna_fm_token

    return df1
    
      
def pickle_save(path, data):
    with open(path, "wb") as f:
        pickle.dump(data, f)

def pickle_read(path):
    with open(path, "rb") as f:
        res = pickle.load(f)      
    return res
    
def pickle_dump(obj, path):
    with open(path, mode='wb') as f:
        pickle.dump(obj,f)

def pickle_load(path):
    with open(path, mode='rb') as f:
        data = pickle.load(f)
        return data    


if __name__ == '__main__':

    start = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--intrain', help='Path')
    parser.add_argument('-it', '--intest', help='Path')
    parser.add_argument('-o', '--outpath', help='Path')
    parser.add_argument('-dm', '--machine', help='Path')
    parser.add_argument('-en', '--encode', help='Path')
    parser.add_argument('-f', '--fold', type=int, help='Path')
    parser.add_argument('-sw', '--seqwin', type=int, help='Path')
    parser.add_argument('-w2v', '--w2vmodel', help='Path')
    
    path = parser.parse_args().intrain
    test_file = parser.parse_args().intest
    out_path_0 = parser.parse_args().outpath
    machine_method = parser.parse_args().machine
    encode_method = parser.parse_args().encode
    kfold = parser.parse_args().fold
    seqwin = parser.parse_args().seqwin
    w2v_path = parser.parse_args().w2vmodel
    
    os.makedirs(out_path_0 + '/' + machine_method + '/' + encode_method, exist_ok=True)
    out_path =  out_path_0 + '/' + machine_method + '/' + encode_method
          
    for i in range(1, kfold+1):
        os.makedirs(out_path + "/" + str(i) + "/data_model", exist_ok=True)
        modelname= "machine_model.sav"
        
        #w2v_model = word2vec.Word2Vec.load(w2v_path + "/" + "rna_w2v.pt")	#W2V
        w2v_model = word2vec.Word2Vec.load(w2v_path + "/" + "rna_w2v_3_100_10_20_1.pt")
        kmer_w2v = 3
        
        train_dataset = pad_input_csv(path + "/" + str(i) + "/cv_train_" + str(i) + ".csv", seqwin, index_col = None) #'seq', 'label'
        val_dataset = pad_input_csv(path + "/" + str(i) + "/cv_val_" + str(i) + ".csv", seqwin, index_col = None)
        test_dataset = pad_input_csv(test_file, seqwin, index_col = None)      

        kw = {}
        train_seq_label=train_dataset.values.tolist()
        val_seq_label  =val_dataset.values.tolist()
        test_seq_label =test_dataset.values.tolist()
        
        if 'W2V' in encode_method:
            train_seq = train_dataset['seq'].tolist()
            valid_seq = val_dataset['seq'].tolist()
            test_seq = test_dataset['seq'].tolist()
            train_X = emb_seq_w2v(train_seq, w2v_model, kmer_w2v)
            valid_X = emb_seq_w2v(valid_seq, w2v_model, kmer_w2v)       
            test_X = emb_seq_w2v(test_seq, w2v_model, kmer_w2v)
            
        elif encode_method == 'Kmer':
            train_X = np.array(Kmer(train_seq_label, k=4, **kw), dtype=float)
            valid_X = np.array(Kmer(val_seq_label, k=4, **kw), dtype=float)
            test_X = np.array(Kmer(test_seq_label, k=4, **kw), dtype=float)
                        
        elif encode_method == 'RCKmer':
            train_X = np.array(RCKmer(train_seq_label, **kw), dtype=float) #k=2
            valid_X = np.array(RCKmer(val_seq_label, **kw), dtype=float)
            test_X = np.array(RCKmer(test_seq_label, **kw), dtype=float)

        elif encode_method == 'DNC':
            train_X = np.array(DNC(train_seq_label, **kw), dtype=float) 
            valid_X = np.array(DNC(val_seq_label, **kw), dtype=float)
            test_X = np.array(DNC(test_seq_label, **kw), dtype=float)               

        elif encode_method == 'TNC':
            train_X = np.array(TNC(train_seq_label, **kw), dtype=float) 
            valid_X = np.array(TNC(val_seq_label, **kw), dtype=float)
            test_X = np.array(TNC(test_seq_label, **kw), dtype=float) 
         
        elif encode_method == 'ENAC':
            kw = {'order': 'ACGU'}
            train_X = np.array(ENAC(train_seq_label, **kw), dtype=float) 
            valid_X = np.array(ENAC(val_seq_label, **kw), dtype=float)
            test_X = np.array(ENAC(test_seq_label, **kw), dtype=float)  

        elif encode_method == 'binary':
            train_X = np.array(binary(train_seq_label, **kw), dtype=float) 
            valid_X = np.array(binary(val_seq_label, **kw), dtype=float)
            test_X = np.array(binary(test_seq_label, **kw), dtype=float)  

        elif encode_method == 'CKSNAP':
            kw = {'order': 'ACGU'}
            train_X = np.array(CKSNAP(train_seq_label, **kw), dtype=float) 
            valid_X = np.array(CKSNAP(val_seq_label, **kw), dtype=float)
            test_X = np.array(CKSNAP(test_seq_label, **kw), dtype=float)  

        elif encode_method == 'NCP':
            train_X = np.array(NCP(train_seq_label, **kw), dtype=float) 
            valid_X = np.array(NCP(val_seq_label, **kw), dtype=float)
            test_X = np.array(NCP(test_seq_label, **kw), dtype=float)  

        elif encode_method == 'ANF':
            train_X = np.array(ANF(train_seq_label, **kw), dtype=float) 
            valid_X = np.array(ANF(val_seq_label, **kw), dtype=float)
            test_X = np.array(ANF(test_seq_label, **kw), dtype=float)  

        elif encode_method == 'PSTNPss': 
            train_X = np.array(PSTNPss(train_seq_label, **kw), dtype=float) 
            valid_X = np.array(PSTNPss(val_seq_label, **kw), dtype=float)
            test_X = np.array(PSTNPss(test_seq_label, **kw), dtype=float)              

        elif encode_method == 'EIIP':
            train_X = np.array(EIIP(train_seq_label, **kw), dtype=float) 
            valid_X = np.array(EIIP(val_seq_label, **kw), dtype=float)
            test_X = np.array(EIIP(test_seq_label, **kw), dtype=float)                                                                           

        elif encode_method == 'PseEIIP':
            train_X = np.array(PseEIIP(train_seq_label, **kw), dtype=float) 
            valid_X = np.array(PseEIIP(val_seq_label, **kw), dtype=float)
            test_X = np.array(PseEIIP(test_seq_label, **kw), dtype=float)                                                                           
        
        elif encode_method == 'CTD':
            train_X = np.array(CTD(train_seq_label, **kw), dtype=float) #k=2
            valid_X = np.array(CTD(val_seq_label, **kw), dtype=float)
            test_X = np.array(CTD(test_seq_label, **kw), dtype=float)
        
        elif encode_method == 'BPB':
            train_X = np.array(BPB(train_seq_label, train_dataset, **kw), dtype=float) 
            valid_X = np.array(BPB(val_seq_label, train_dataset, **kw), dtype=float)
            test_X = np.array(BPB(test_seq_label, train_dataset, **kw), dtype=float)               

        elif encode_method == 'NCP_ND':
            train_X = np.array(NCP_ND(train_seq_label, **kw), dtype=float) 
            valid_X = np.array(NCP_ND(val_seq_label, **kw), dtype=float)
            test_X = np.array(NCP_ND(test_seq_label, **kw), dtype=float) 
                         
        elif encode_method == 'NPS':
            kw = {'order': 'ACGU'}
            train_X = np.array(NPS(train_seq_label, **kw), dtype=float) 
            valid_X = np.array(NPS(val_seq_label, **kw), dtype=float)
            test_X = np.array(NPS(test_seq_label, **kw), dtype=float)  
                    
        elif encode_method == 'NPPS':
            train_X = np.array(NPPS(train_seq_label, train_dataset, xi=0, **kw), dtype=float) 
            valid_X = np.array(NPPS(val_seq_label, train_dataset, xi=0, **kw), dtype=float)
            test_X = np.array(NPPS(test_seq_label, train_dataset, xi=0, **kw), dtype=float)  

        elif encode_method == 'PseKNC':
            kw = {'order': 'ACGU'}
            train_X = np.array(PseKNC(train_seq_label, k=3, **kw), dtype=float) 
            valid_X = np.array(PseKNC(val_seq_label, k=3, **kw), dtype=float)
            test_X = np.array(PseKNC(test_seq_label, k=3, **kw), dtype=float)
                                                                                                                                                                                                                  
        else :
            pass
            print('no encode method')
            exit()
          
         
        train_y = train_dataset['label'].to_numpy()    
        valid_y = val_dataset['label'].to_numpy() 
        test_y = test_dataset['label'].to_numpy()
            
        cv_result = np.zeros((len(valid_y), 2))
        cv_result[:, 1] = valid_y
        test_result = np.zeros((len(test_y), 2))
        test_result[:,1] = test_y   #score:one of two, label
     
        if machine_method == 'RF':
            model = RandomForestClassifier(max_depth=4, random_state=0, n_estimators=100)
            clf = model.fit(train_X, train_y)

        elif machine_method == 'NB':
            model = GaussianNB() # 正規分布を仮定したベイズ分類
            clf = model.fit(train_X, train_y)
                    
        elif machine_method == 'KNN':
            model = KNeighborsClassifier()
            clf = model.fit(train_X, train_y)

        elif machine_method == 'LR':
            model = LogisticRegression(random_state=0)
            clf = model.fit(train_X, train_y)  
                    
        elif machine_method == 'SVM':    
            model = svm.SVC(probability=True)
            clf = model.fit(train_X, train_y)

        elif machine_method == 'XGB':
            xgb_train = xgb.DMatrix(train_X, train_y)
            xgb_eval  = xgb.DMatrix(valid_X , valid_y)
            params = {
                "learning_rate": 0.01,
                "max_depth": 3
            }
            clf = xgb.train(params, 
              xgb_train, evals=[(xgb_train, "train"), (xgb_eval, "validation")], 
              num_boost_round=100, early_stopping_rounds=20)
        
        elif machine_method == 'CBC':
            model = CatBoostClassifier(eval_metric = 'AUC', loss_function = "Logloss", task_type="CPU", learning_rate=0.01, depth=9, iterations=500)
            clf = model.fit(train_X, train_y) 
            
        elif machine_method == 'LGBM':   
            lgb_train = lgb.Dataset(train_X , train_y)
            lgb_eval = lgb.Dataset(valid_X , valid_y, reference=lgb_train)
            params = {         
                'objective': 'binary',        
                'metric': 'auc',         
                'verbosity': -1,
                'random_state': 123,
            }
            clf = lgb.train(params, lgb_train, valid_sets=lgb_eval,
                      verbose_eval=50,
                      num_boost_round=1000,
                      early_stopping_rounds=100) 
        else:
            print('No learning method')
            exit()

        pickle.dump(clf, open(out_path + "/" + str(i) + "/data_model/machine_model.asv",'wb'))  
     
        #CV
        if machine_method == 'LGBM':  
            score = clf.predict(valid_X, num_iteration=clf.best_iteration)
            cv_result[:, 0] = score
        elif machine_method == 'XGB':  
            score = clf.predict(xgb_eval)
            cv_result[:, 0] = score
        else :
            score = clf.predict_proba(valid_X)
            cv_result[:, 0] = score[:,1]

        #independent test
        if test_dataset.shape[0] != 0:
            if machine_method == 'LGBM':
                test_result[:, 0] = clf.predict(test_X, num_iteration=clf.best_iteration)
            elif machine_method == 'XGB': 
                test_result[:, 0] = clf.predict(xgb.DMatrix(test_X))
            else:
                test_result[:, 0] = clf.predict_proba(test_X)[:,1]
      
        #CV  
        cv_output = pd.DataFrame(cv_result,  columns=['prob', 'label'] )
        cv_output.to_csv(out_path  + "/" + str(i) + "/val_roc.csv")  #prob, label

        #independent test
        test_output = pd.DataFrame(test_result,  columns=['prob', 'label'] )
        test_output.to_csv(out_path  + "/" + str(i) + "/test_roc.csv")  #prob, label

    print('elapsed time', time.time() - start)

