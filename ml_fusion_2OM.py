import os
import sys
import pickle
import pandas as pd
import numpy as np
import time
import argparse
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
sys.path.append(os.path.abspath(''))
sys.path.append(os.path.abspath('..'))
from valid_metrices_p2 import *

item_column =['Thre','Rec','Pre','F1', 'Spe', 'Acc','MCC','AUC','PRAUC'] 

def combine_model(valid_data, test_data, data_path, out_dir, kfold):  
    prediction_result_cv = []
    prediction_result_test = []
    
    train_y, train_X = valid_data[:,0], valid_data[:,1:]
    valid_y, valid_X = valid_data[:,0], valid_data[:,1:]
    test_y, test_X = test_data[:,0], test_data[:,1:]
    
    model = LogisticRegression(random_state=0)
    rfc = model.fit(train_X, train_y)
    
    if os.path.isfile('%s/result/combine/%s/%s/lr_model.asv' % (data_path, out_dir, kfold)) :
        pass
    else:
       os.makedirs('%s/result/combine/%s/%s' % (data_path, out_dir, kfold), exist_ok=True)
    pickle.dump(rfc, open('%s/result/combine/%s/%s/lr_model.asv' % (data_path, out_dir, kfold), 'wb'))

    
    """
    if os.path.isfile('%s/result/%s/%s/lr_model.asv' % (data_path, out_dir, kfold)):
        rfc = pickle.load(open('%s/result/%s/%s/lr_model.asv' % (data_path, out_dir, kfold), 'rb'))
    else :
        print('No %s model' %out_dir)
    """
                 
    scores = rfc.predict_proba(valid_X)
    tmp_result = np.zeros((len(valid_y), 1+2))
    tmp_result[:, 0], tmp_result[:, 1:] = valid_y, scores    
    prediction_result_cv.append(tmp_result)
             
    if test_X.shape[0] != 0:
        scores_test = rfc.predict_proba(test_X)
        tmp_result_test = np.zeros((len(test_y), 1+2))
        tmp_result_test[:, 0], tmp_result_test[:, 1:] = test_y, scores_test    
        prediction_result_test.append(tmp_result_test)           
        #print(prediction_result_cv[0]) #[1  0.094, 0.9057 ]
    
    valid_probs = prediction_result_cv[0][:,2]
    valid_labels = prediction_result_cv[0][:,0]    
    #print(np.array([valid_probs, valid_labels]).T)
    
    test_probs = prediction_result_test[0][:,2]
    test_labels = prediction_result_test[0][:,0]   

    cv_output = pd.DataFrame(np.array([valid_probs, valid_labels]).T,  columns=['prob', 'label'] )
    cv_output.to_csv('%s/result/combine/%s/%s/val_roc.csv' % (data_path, out_dir, kfold))   
    test_output = pd.DataFrame(np.array([test_probs, test_labels]).T,  columns=['prob', 'label'] )
    test_output.to_csv('%s/result/combine/%s/%s/test_roc.csv' % (data_path, out_dir, kfold))
    
    #print(f'validation: prob label {valid_probs} {valid_labels} ')
    
    # metrics calculation
    th_, rec_, pre_, f1_, spe_, acc_, mcc_, auc_, pred_class, prauc_ = eval_metrics(valid_probs, valid_labels) 
    valid_matrices = th_, rec_, pre_, f1_, spe_, acc_, mcc_, auc_, prauc_
    th_, rec_, pre_, f1_, spe_, acc_, mcc_, auc_, pred_class, prauc_ = th_eval_metrics(th_, test_probs, test_labels)
    test_matrices = th_, rec_, pre_, f1_, spe_, acc_, mcc_, auc_, prauc_

    print_results(valid_matrices, test_matrices) 
    
    #print(f'valid_matrices {valid_matrices}')  
    
    
    df = pd.DataFrame([valid_matrices, test_matrices], index=['valid','test'], columns=item_column)
    df2 = pd.DataFrame([test_matrices], index=['test'], columns=item_column)
    
    return df, df2


def train_test(kfold, data_path, out_dir, combination):
    #feature combine for each fold
    valid_data =[]
    test_data =[]
    for comb in combination:
        machine = comb[0]
        fea = comb[1]
        for datype in ['val','test']:
                fea_file = data_path + '/result/%s/%s/%s/%s_roc.csv' %(machine, fea, str(kfold), datype)
                fea_data = pd.read_csv(fea_file)
                if datype =='val':
                    valid_data.append(fea_data['label'].values.tolist())
                    valid_data.append(fea_data['prob'].values.tolist())
                elif datype =='test':
                    test_data.append(fea_data['label'].values.tolist())
                    test_data.append(fea_data['prob'].values.tolist())

    valid_data = np.array(valid_data).T
    test_data = np.array(test_data).T    
    #print(valid_data) 
    print('valid_data.shape {}'.format(valid_data.shape) )
    print(f'valid_data[:,0]: {valid_data[:,0]}')
    print(f'valid_data[:,1:]: {valid_data[:,1:]}')
    
    # Redundant labels [label,prob,label, prob,....] are removed 
    valid_data = np.delete(valid_data, [i for i in range(2, 2*len(combination), 2)], 1)
    test_data  = np.delete(test_data, [i for i in  range(2, 2*len(combination), 2)], 1)           
    print('non redundant valid_data.shape {}'.format(valid_data.shape) )
    
    # training and testing
    df, df2 = combine_model(valid_data, test_data, data_path, out_dir, kfold)

    return df, df2
    

# score combine method based on logistic regression
if __name__ == '__main__':
    
    kfold=5
    out_name ='combine'
    outfile = out_name +'.csv'

    os.chdir("..")
    os.chdir("..")
    source_path = os.getcwd()
    data_path = source_path +"/data"
    print(data_path)

    machine_method = ["LGBM","CBC","XGB","RF","SVM","LR","NB","KNN"]
    encode_method = ["binary","NCP","EIIP","Kmer", "RCKmer","DNC","TNC","CKSNAP","ANF","PseEIIP","ENAC","CTD","BPB","NCP_ND","NPS","NPPS","PseKNC","W2V"]
          
    combination = []
    for machine in machine_method:
        for encode in encode_method:
            combination.append([machine, encode])
        
       
    outfile = data_path +'/result/%s' %outfile

    for k in range(1, kfold+1):
        df, df2 = train_test(k, data_path, out_name, combination)
        if k==1:        
            df_cat = df    #train
            df_cat_2 = df2 #test
            #print(df2)
        else :
            df_cat = pd.concat([df_cat, df], axis=0)
            df_cat_2 = pd.concat([df_cat_2, df2], axis=0)

    
    #df_cat.to_csv(outfile)
    print(df_cat.mean())
    print(df_cat_2.mean()) 
    
    df_cat_21 = pd.DataFrame(columns=item_column)
    df_cat_21.loc["mean_train"] = df_cat.mean()
    df_cat_21.loc["sd_train"] = df_cat.std()
    df_cat_21.loc["mean_test"] = df_cat_2.mean()
    df_cat_21.loc["sd_test"] = df_cat_2.std()
    
    df_cat_21.to_csv(outfile)
    print(df_cat_21 )
    



    
    
    
       
