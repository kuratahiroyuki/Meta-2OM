#!/bin/bash

current_path=$(pwd)
echo ${current_path}
cd ..
cd ..
main_path=$(pwd)
cd ${current_path}
echo $(pwd)


train_path=${main_path}/data/dataset/cross_val
test_file=${main_path}/data/dataset/independent_test/independent_test.csv
result_path=${main_path}/data/result
w2v_path=${main_path}/data/w2v_model

kfold=5
seqwin=41

for machine_method in LGBM CBC XGB RF SVM LR NB KNN
do
for encode_method in binary NCP EIIP Kmer RCKmer DNC TNC CKSNAP ANF PseEIIP ENAC CTD BPB NCP_ND NPS NPPS PseKNC W2V
do 
echo ${machine_method}
echo ${encode_method}

python ml_train_test_2OM.py  --intrain ${train_path} --intest ${test_file} --outpath ${result_path} --machine ${machine_method}  --encode ${encode_method} --fold ${kfold} --seqwin ${seqwin} --w2vmodel ${w2v_path}

done
done

