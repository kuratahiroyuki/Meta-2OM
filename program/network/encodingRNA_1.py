#!/usr/bin/env python
#_*_coding:utf-8_*_

import pandas as pd
import numpy as np
from collections import Counter
import re, sys, os, platform
import itertools
from kmer import get_kmers, kmer2number, get_dinucNum_with_interval, number2kmer
from get_RNAPhyche import get_RNAPhyche 
#from check_sequences import *

###########################################33
def check_fasta_with_equal_length(fastas):
    status = True
    lenList = set()
    for i in fastas:
        lenList.add(len(i[0])) ### i[1] -> i[0] に変更
    if len(lenList) == 1:
        return True
    else:
        return False

def get_min_sequence_length(fastas):
    minLen = 10000
    for i in fastas:
        if minLen > len(i[0]): ###
            minLen = len(i[0]) ###
    return minLen

def get_min_sequence_length_1(fastas):
    minLen = 10000
    for i in fastas:
        if minLen > len(re.sub('-', '', i[0])): ###
            minLen = len(re.sub('-', '', i[0])) ###
    return minLen


def get_min_sequence_length(fastas):
    minLen = 10000
    for i in fastas:
        if minLen > len(i[0]): ###
            minLen = len(i[0]) ###
    return minLen

def get_min_sequence_length_1(fastas):
    minLen = 10000
    for i in fastas:
        if minLen > len(re.sub('-', '', i[0])): ###
            minLen = len(re.sub('-', '', i[0])) ###
    return minLen

###
def kmerArray(sequence, k):
    kmer = []
    for i in range(len(sequence) - k + 1):
        kmer.append(sequence[i:i + k])
    return kmer

nuc2num_dict = {'A': 0, 'C': 1, 'G': 2, 'U': 3}
num2nuc_dict = {0: 'A', 1: 'C', 2: 'G', 3: 'U'}

##############################################################################
def EIIP(fastas, **kw):
    if check_fasta_with_equal_length == False:
        print('Error: for "EIIP" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    AA = 'ACGU'

    EIIP_dict ={
        'A': 0.1260,
        'C': 0.1340,
        'G': 0.0806,
        'U': 0.1335,
        '-': 0,
    }

    encodings = []
    header = ['#', 'label']

    for i in range(1, len(fastas[0][0]) + 1): ### fastas[0][1] => fastas[0][0]
        header.append('F'+str(i))
    #encodings.append(header)

    for i in fastas:
        #name, sequence, label = i[0], i[1], i[2]
        #code = [name, label]
        sequence, label = i[0], i[1]
        code = []
        for aa in sequence:
            code.append(EIIP_dict.get(aa, 0))
        encodings.append(code)
    return encodings


def TriNcleotideComposition(sequence, base):
    trincleotides = [nn1 + nn2 + nn3 for nn1 in base for nn2 in base for nn3 in base]
    tnc_dict = {}
    for triN in trincleotides:
        tnc_dict[triN] = 0
    for i in range(len(sequence) - 2):
        tnc_dict[sequence[i:i + 3]] += 1
    for key in tnc_dict:
       tnc_dict[key] /= (len(sequence) - 2)
    return tnc_dict


def PseEIIP(fastas, **kw):
    for i in fastas:
        if re.search('[^ACGU-]', i[0]):###
            print('Error: illegal character included in the fasta sequences, only the "ACGU-" are allowed by this PseEIIP scheme.')
            return 0

    base = 'ACGU'

    EIIP_dict = {
        'A': 0.1260,
        'C': 0.1340,
        'G': 0.0806,
        'U': 0.1335,
    }

    trincleotides = [nn1 + nn2 + nn3 for nn1 in base for nn2 in base for nn3 in base]
    EIIPxyz = {}
    for triN in trincleotides:
        EIIPxyz[triN] = EIIP_dict[triN[0]] + EIIP_dict[triN[1]] + EIIP_dict[triN[2]]

    encodings = []
    header = ['#', 'label'] + trincleotides
    #encodings.append(header)

    for i in fastas:
        #name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        sequence, label =  re.sub('-', '', i[0]), i[1]
        #code = [name, label]
        code = []
        trincleotide_frequency = TriNcleotideComposition(sequence, base)
        code = code + [EIIPxyz[triN] * trincleotide_frequency[triN] for triN in trincleotides]
        encodings.append(code)
    return encodings


##############################################################################
def CTD(fastas, **kw):
    
    base = 'ACGU'
    encodings = []
    
    for i in fastas:
        #name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        seq, label = i[0], i[1]
        #code = [name, label]
                
        n = float(len(seq) - 1)
        num_A, num_U, num_G, num_C = 0.0, 0.0, 0.0, 0.0
        AU_trans, AG_trans, AC_trans, UG_trans, UC_trans, GC_trans = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        for i in range(len(seq) - 1):
            if seq[i] == "A":
                num_A = num_A + 1
            if seq[i] == "U":
                num_U = num_U + 1
            if seq[i] == "G":
                num_G = num_G + 1
            if seq[i] == "C":
                num_C = num_C + 1
            if (seq[i] == "A" and seq[i + 1] == "U") or (seq[i] == "U" and seq[i + 1] == "A"):
                AU_trans = AU_trans + 1
            if (seq[i] == "A" and seq[i + 1] == "G") or (seq[i] == "G" and seq[i + 1] == "A"):
                AG_trans = AG_trans + 1
            if (seq[i] == "A" and seq[i + 1] == "C") or (seq[i] == "C" and seq[i + 1] == "A"):
                AC_trans = AC_trans + 1
            if (seq[i] == "U" and seq[i + 1] == "G") or (seq[i] == "G" and seq[i + 1] == "U"):
                UG_trans = UG_trans + 1
            if (seq[i] == "U" and seq[i + 1] == "C") or (seq[i] == "C" and seq[i + 1] == "U"):
                UC_trans = UC_trans + 1
            if (seq[i] == "G" and seq[i + 1] == "C") or (seq[i] == "C" and seq[i + 1] == "G"):
                GC_trans = GC_trans + 1

        a, u, g, c = 0, 0, 0, 0
        A0_dis, A1_dis, A2_dis, A3_dis, A4_dis = 0.0, 0.0, 0.0, 0.0, 0.0
        U0_dis, U1_dis, U2_dis, U3_dis, U4_dis = 0.0, 0.0, 0.0, 0.0, 0.0
        G0_dis, G1_dis, G2_dis, G3_dis, G4_dis = 0.0, 0.0, 0.0, 0.0, 0.0
        C0_dis, C1_dis, C2_dis, C3_dis, C4_dis = 0.0, 0.0, 0.0, 0.0, 0.0
        for i in range(len(seq) - 1):
            if seq[i] == "A":
                a = a + 1
                if a == 1:
                    A0_dis = ((i * 1.0) + 1) / n
                if a == int(round(num_A / 4.0)):
                    A1_dis = ((i * 1.0) + 1) / n
                if a == int(round(num_A / 2.0)):
                    A2_dis = ((i * 1.0) + 1) / n
                if a == int(round((num_A * 3 / 4.0))):
                    A3_dis = ((i * 1.0) + 1) / n
                if a == num_A:
                    A4_dis = ((i * 1.0) + 1) / n
            if seq[i] == "U":
                u = u + 1
                if u == 1:
                    U0_dis = ((i * 1.0) + 1) / n
                if u == int(round(num_U / 4.0)):
                    U1_dis = ((i * 1.0) + 1) / n
                if u == int(round((num_U / 2.0))):
                    U2_dis = ((i * 1.0) + 1) / n
                if u == int(round((num_U * 3 / 4.0))):
                    U3_dis = ((i * 1.0) + 1) / n
                if u == num_U:
                    U4_dis = ((i * 1.0) + 1) / n
            if seq[i] == "G":
                g = g + 1
                if g == 1:
                    G0_dis = ((i * 1.0) + 1) / n
                if g == int(round(num_G / 4.0)):
                    G1_dis = ((i * 1.0) + 1) / n
                if g == int(round(num_G / 2.0)):
                    G2_dis = ((i * 1.0) + 1) / n
                if g == int(round(num_G * 3 / 4.0)):
                    G3_dis = ((i * 1.0) + 1) / n
                if g == num_G:
                    G4_dis = ((i * 1.0) + 1) / n
            if seq[i] == "C":
                c = c + 1
                if c == 1:
                    C0_dis = ((i * 1.0) + 1) / n
                if c == int(round(num_C / 4.0)):
                    C1_dis = ((i * 1.0) + 1) / n
                if c == int(round(num_C / 2.0)):
                    C2_dis = ((i * 1.0) + 1) / n
                if c == int(round(num_C * 3 / 4.0)):
                    C3_dis = ((i * 1.0) + 1) / n
                if c == num_C:
                    C4_dis = ((i * 1.0) + 1) / n
        CTD_feature = [num_A / n, num_U / n, num_G / n, num_C / n,
                       AU_trans / n - 1, AG_trans / (n - 1), AC_trans / (n - 1),
                       UG_trans / n - 1, UC_trans / (n - 1), GC_trans / (n - 1),
                       A0_dis, A1_dis, A2_dis, A3_dis, A4_dis,
                       U0_dis, U1_dis, U2_dis, U3_dis, U4_dis,
                       G0_dis, G1_dis, G2_dis, G3_dis, G4_dis,
                       C0_dis, C1_dis, C2_dis, C3_dis, C4_dis]
        
        encodings.append(CTD_feature)
    
    return encodings

##############################################################################
d_max = 3

bn_dict = {('A', 'A'): 0, ('A', 'U'): 1, ('A', 'C'): 2, ('A', 'G'): 3, ('U', 'A'): 4, ('U', 'U'): 5, ('U', 'C'): 6,
           ('U', 'G'): 7, ('C', 'A'): 8, ('C', 'U'): 9, ('C', 'C'): 10, ('C', 'G'): 11, ('G', 'A'): 12, ('G', 'U'): 13,
           ('G', 'C'): 14, ('G', 'G'): 15}

def NPS(fastas, **kw):
    
    base = 'ACGU'
    encodings = []
    
    for i in fastas:
        #name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        seq, label = i[0], i[1]
        #code = [name, label]
        res = []
        seq_len = len(seq)
        if seq_len < d_max + 2:
            print("error")
            return None
        for d in range(d_max + 1):
            pair_count = np.zeros(16)
            for i in range(seq_len - 1 - d):
                index = bn_dict[(seq[i], seq[i + d + 1])]
                pair_count[index] += 1
            res.extend(pair_count / (seq_len - d - 1))
        
        encodings.append(res)
    
    return encodings

##############################################################################
cpm_dict = {'A': (1., 1., 1.), 'C': (0., 0., 1.), 'G': (1., 0., 0.), 'U': (0., 1., 0.)}
cpm_dict2 = {'A': [0., 0., 0., 1.], 'C': [0., 0., 1., 0.], 'G': [0., 1., 0., 0.], 'U': [1., 0., 0., 0.]}
seq_dict = {'A': 0, 'C': 1, 'G': 2, 'U': 3}


def density(seq):
    res = []
    d = {'A': 0., 'U': 0., 'C': 0., 'G': 0.}
    for i in range(len(seq)):
        d[seq[i]] += 1
        res.append(d[seq[i]] / (i + 1))
    return res


def NCP_ND(fastas, **kw):
    
    base = 'ACGU'
    encodings = []
    
    for i in fastas:
        #name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        seq, label = i[0], i[1]
        #code = [name, label]
        res = []
        # seq = remove_center_GAC(seq)
        den = density(seq)
        for n, i in zip(seq, range(len(den))):
            res.extend(cpm_dict[n])
            res.append(den[i])
        
        encodings.append(res)
    
    return encodings

##############################################################################
def get_pos_neg_subsets_df(data_df):  # Gets a positive samples subset and negative samples subset in the DataFrame
    #data_df = pd.DataFrame(train_seq_label, columns=["seq","label"])
    positive_subsets_df = data_df[data_df.label == 1]
    negative_subsets_df = data_df[data_df.label == 0]
    
    return positive_subsets_df, negative_subsets_df

def get_kmer_position_frequency(data_df, k):
    seq_df = data_df['seq'].apply(get_kmers, k=k).apply(pd.Series)
    position_frequency = np.zeros((4 ** k, seq_df.columns.size))
    for col in seq_df.columns:
        # print(col)
        kmer_p = seq_df[col].value_counts(normalize=True)  # 某个位置的kmer的出现频率
        for kmerN, p in zip(kmer_p.index, kmer_p):
            # print(nuc, p)
            position_frequency[int(kmerN)][col] = p
    return position_frequency

def get_pos_neg_posteriori_probability(data_df):  
    positive_subsets_df, negative_subsets_df = get_pos_neg_subsets_df(data_df)

    positive_posteriori_probability = get_kmer_position_frequency(positive_subsets_df, k=1)
    negative_posteriori_probability = get_kmer_position_frequency(negative_subsets_df, k=1)
    return positive_posteriori_probability, negative_posteriori_probability

def BPB(fastas, train_data_df, **kw):
    
    base = 'ACGU'
    encodings = []
    
    positive_posteriori_probability, negative_posteriori_probability = get_pos_neg_posteriori_probability(train_data_df)
        
    for i in fastas:
        #name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        seq, label = i[0], i[1]
        #code = [name, label]
        BPB_feature = []
        if 'N' in seq:
            print("error! The N can't in seq.")
        if len(seq) != positive_posteriori_probability.shape[1] or len(seq) != negative_posteriori_probability.shape[1]:
            print(
                "error! The length of seq is not equal to positive_posteriori_probability or negative_posteriori_probability.")
        for i in range(len(seq)):
            BPB_feature.append(positive_posteriori_probability[nuc2num_dict[seq[i]]][i])
            BPB_feature.append(negative_posteriori_probability[nuc2num_dict[seq[i]]][i])
        
        encodings.append(BPB_feature)
    
    return encodings

##############################################################################
def get_dinuc_with_interval_postion_frequency(data_df, xi):
    seq_df = data_df['seq'].apply(get_dinucNum_with_interval, xi=xi).apply(pd.Series)
    position_frequency = np.zeros((4 ** 2, seq_df.columns.size))
    for col in seq_df.columns:
        kmer_p = seq_df[col].value_counts(normalize=True)  # 某个位置的kmer出现的频率
        for kmerN, p in zip(kmer_p.index, kmer_p):
            # print(nuc, p)
            # print(kmerN, p, col)
            position_frequency[int(kmerN)][col] = p
    return position_frequency

def get_front_post_Tp(data_df, xi):
    Ts = get_kmer_position_frequency(data_df, k=1)
    Td = get_dinuc_with_interval_postion_frequency(data_df, xi=xi)
    front_Tp = np.zeros(Td.shape)
    post_Tp = np.zeros(Td.shape)
    for i in range(front_Tp.shape[0]):
        for j in range(front_Tp.shape[1]):
            if Ts[int(i / 4)][j] == 0:
                front_Tp[i][j] = 0.
            else:
                front_Tp[i][j] = Td[i][j] / Ts[int(i / 4)][j]
            if Ts[int(i / 4)][j + xi + 1] == 0:
                post_Tp[i][j] = 0.
            else:
                post_Tp[i][j] = Td[i][j] / Ts[int(i / 4)][j + xi + 1]
    return front_Tp, post_Tp

def get_pos_neg_front_post_Tp(data_df, xi):
    positive_subsets_df, negative_subsets_df = get_pos_neg_subsets_df(data_df)
    pos_front_Tp, pos_post_Tp = get_front_post_Tp(positive_subsets_df, xi=xi)
    neg_front_Tp, neg_post_Tp = get_front_post_Tp(negative_subsets_df, xi=xi)
    return pos_front_Tp, pos_post_Tp, neg_front_Tp, neg_post_Tp


def NPPS(fastas, train_data_df, xi=3, **kw):
    
    base = 'ACGU'
    encodings = []
    
    pos_front_Tp, pos_post_Tp, neg_front_Tp, neg_post_Tp = get_pos_neg_front_post_Tp(train_data_df, xi)
    
    
    for i in fastas:
        #name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        seq, label = i[0], i[1]
        #code = [name, label]
        NPPS_feature = []
        dinucNums = get_dinucNum_with_interval(seq, xi=xi)
        for j, dinucNum in zip(range(len(dinucNums)), dinucNums):
            # NPPS_feature.append(pos_front_Tp[dinucNum][i])
            # NPPS_feature.append(pos_post_Tp[dinucNum][i])
            # NPPS_feature.append(neg_front_Tp[dinucNum][i])
            # NPPS_feature.append(neg_post_Tp[dinucNum][i])
            NPPS_feature.append(pos_post_Tp[dinucNum][j] - neg_post_Tp[dinucNum][j])
        # print(NPPS_feature)
        encodings.append(NPPS_feature)
    
    return encodings

########################################################################################

def get_NC_frequency_from_seq(seq, k, type='list'):
    '''
    The frequency of kmers is counted from the sequence.
    :param seq: input sequence
    :param k:
    :param type: if the value is "list"，then return a list of frequency according to the kmerNumber; else return a dict.
    :return:
    '''
    if type != 'list' and type != 'dict':
        print("error! Return normalize_type must be list or dict")

    counter = Counter(get_kmers(seq, k))
    frequence_list = [0 for i in range(4 ** k)]
    for key in counter.keys():
        frequence_list[key] += counter[key]
    frequence_list = list(np.array(frequence_list) / np.array(frequence_list).sum())
    if type == 'list':
        return frequence_list
    else:
        frequence_dict = dict()
        for i in range(4 ** k):
            frequence_dict[number2kmer(i, k=k)] = frequence_list[i]
        return frequence_dict

def get_correlationValue(diPC_df):
    correlationValue = np.zeros((len(diPC_df.columns), len(diPC_df), len(diPC_df)))
    for i in range(len(diPC_df.columns)):
        for j in range(len(diPC_df)):
            for k in range(len(diPC_df)):
                correlationValue[i, j, k] = (diPC_df.iloc[j, i] - diPC_df.iloc[k, i]) ** 2
    return correlationValue


def get_theta_array(seq, lambde, correlationValue):
    theta_array = []
    kmers = get_kmers(seq, 2)
    for ilambda in range(1, lambde + 1):
        theta = 0
        # for i in range(len(seq) - ilambda - 1):
        for i in range(len(seq) - lambde - 1):  # repDNA做法
            pepA = kmers[i]
            pepB = kmers[i + ilambda]

            CC = 0
            for j in range(len(correlationValue)):
                # print(j, pepA, pepB)
                CC += correlationValue[j, pepA, pepB]
                # print(correlationValue[j, pepA, pepB])
            CC /= len(correlationValue)
            theta += CC
        theta_array.append(theta / (len(seq) - ilambda - 1))
    return theta_array

def get_PseDNC_feature(seq, weight, lambde, correlationValue):
    '''
    :param seq:
    :param weight:
    :param lambde:
    :return:
    '''
    PseDNC_feature = []
    DNC_frequency = get_NC_frequency_from_seq(seq, k=2)
    theta_array = get_theta_array(seq, lambde, correlationValue=correlationValue)
    for i in range(4 ** 2):
        PseDNC_feature.append(DNC_frequency[i] / (1 + weight * sum(theta_array)))
    for i in range(4 ** 2 + 1, 4 ** 2 + lambde + 1):
        PseDNC_feature.append((weight * theta_array[i - 17]) / (1 + weight * sum(theta_array)))
    return PseDNC_feature
    

def PseKNC(fastas, k=3, weight=0.5, lambde=3, **kw):
    
    base = 'ACGU'
    encodings = []
    
    di_phy_list = ['Rise(RNA)', 'Roll(RNA)', 'Shift(RNA)', 'Slide(RNA)', 'Tilt(RNA)', 'Twist(RNA)']
    diPC_df = get_RNAPhyche(phy_list=di_phy_list, k=2, standardized=False)
    diPC_df['kmer'] = diPC_df.index
    diPC_df['kmer'] = diPC_df['kmer'].apply(kmer2number)
    diPC_df = diPC_df.set_index('kmer')
    # print(diPC_df)
    
    correlationValue = get_correlationValue(diPC_df)
    
    for i in fastas:
        #name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        seq, label = i[0], i[1]
        #code = [name, label]
        PseKNC_feature = []
        KNC_frequency = get_NC_frequency_from_seq(seq, k=k)
        theta_array = get_theta_array(seq, lambde, correlationValue=correlationValue)
        # print(theta_array)
        for i in range(4 ** k):
            PseKNC_feature.append(KNC_frequency[i] / (1 + weight * sum(theta_array)))
        for i in range(4 ** k + 1, 4 ** k + lambde + 1):
            PseKNC_feature.append((weight * theta_array[i - 4 ** k - 1]) / (1 + weight * sum(theta_array)))
        
        encodings.append(PseKNC_feature)
    
    return encodings

#################################################################################################################































