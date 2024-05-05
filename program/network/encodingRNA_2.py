#!/usr/bin/env python
#_*_coding:utf-8_*_

import pandas as pd
import numpy as np
from collections import Counter
import re, sys, os, platform
import itertools
import check_sequences


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

###
def Kmer(fastas, k=2, type="RNA", upto=False, normalize=True, **kw):
    encoding = []
    header = ['#', 'label']
    NA = 'ACGU'
    if type in ("DNA", 'RNA'):
        NA = 'ACGU'
    else:
        NA = 'ACDEFGHIKLMNPQRSTVWY'

    if k < 1:
        print('Error: the k-mer value should larger than 0.')
        return 0

    if upto == True:
        for tmpK in range(1, k + 1):
            for kmer in itertools.product(NA, repeat=tmpK):
                header.append(''.join(kmer))
        #encoding.append(header)
        #print(header)

        for i in fastas:
            #name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
            sequence, label = i[0], i[1]
            count = Counter()
            for tmpK in range(1, k + 1):
                kmers = kmerArray(sequence, tmpK)
                count.update(kmers)
                if normalize == True:
                    for key in count:
                        if len(key) == tmpK:
                            count[key] = count[key] / len(kmers)
            #code = [name, label]
            code = []
            for j in range(2, len(header)):
                if header[j] in count:
                    code.append(count[header[j]])
                else:
                    code.append(0)
            encoding.append(code)
    else:
        for kmer in itertools.product(NA, repeat=k):
            header.append(''.join(kmer))
        #encoding.append(header) #['AA', 'AC', 'AG', 'AU', 'CA', 'CC', 'CG', 'CU', 'GA', 'GC', 'GG', 'GU', 'UA', 'UC', 'UG', 'UU'] #16

        for i in fastas:
            #name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
            sequence, label = i[0], i[1]
            kmers = kmerArray(sequence, k)
            #print(f'kmers {len(kmers)}, {kmers}')
            count = Counter()
            count.update(kmers)
            #print(f'count {count}')
            if normalize == True:
                for key in count:
                    count[key] = count[key] / len(kmers)
            #code = [name, label]
            code = []
            for j in range(2, len(header)):
                if header[j] in count:
                    code.append(count[header[j]])
                else:
                    code.append(0)
      
            encoding.append(code)
            #print(f'code {code}')


    return encoding


def RC(kmer):
    myDict = {
        'A': 'U',
        'C': 'G',
        'G': 'C',
        'U': 'A'
    }
    return ''.join([myDict[nc] for nc in kmer[::-1]])


def generateRCKmer(kmerList):
    rckmerList = set()
    myDict = {
        'A': 'U',
        'C': 'G',
        'G': 'C',
        'U': 'A'
    }
    for kmer in kmerList:
        rckmerList.add(sorted([kmer, ''.join([myDict[nc] for nc in kmer[::-1]])])[0])
    return sorted(rckmerList)

###
def RCKmer(fastas, k=2, upto=False, normalize=True, **kw):
    encoding = []
    header = ['#', 'label']
    NA = 'ACGU'

    if k < 1:
        print('Error: the k-mer value should larger than 0.')
        return 0

    if upto == True:
        for tmpK in range(1, k + 1):
            tmpHeader = []
            for kmer in itertools.product(NA, repeat=tmpK):
                tmpHeader.append(''.join(kmer))
            header = header + generateRCKmer(tmpHeader)
        myDict = {}
        for kmer in header[2:]:
            rckmer = RC(kmer)
            if kmer != rckmer:
                myDict[rckmer] = kmer
        #encoding.append(header)
        for i in fastas:
            #name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
            sequence, label = i[0], i[1]
            count = Counter()
            for tmpK in range(1, k + 1):
                kmers = kmerArray(sequence, tmpK)
                for j in range(len(kmers)):
                    if kmers[j] in myDict:
                        kmers[j] = myDict[kmers[j]]
                count.update(kmers)
                if normalize == True:
                    for key in count:
                        if len(key) == tmpK:
                            count[key] = count[key] / len(kmers)
            #code = [name, label]
            code = []
            for j in range(2, len(header)):
                if header[j] in count:
                    code.append(count[header[j]])
                else:
                    code.append(0)
            encoding.append(code)
    else:
        tmpHeader = []
        for kmer in itertools.product(NA, repeat=k):
            tmpHeader.append(''.join(kmer))
        header = header + generateRCKmer(tmpHeader)
        myDict = {}
        for kmer in header[2:]:
            rckmer = RC(kmer)
            if kmer != rckmer:
                myDict[rckmer] = kmer

        #encoding.append(header)
        for i in fastas:
            #name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
            sequence, label = i[0], i[1]
            kmers = kmerArray(sequence, k)
            for j in range(len(kmers)):
                if kmers[j] in myDict:
                    kmers[j] = myDict[kmers[j]]
            count = Counter()
            count.update(kmers)
            if normalize == True:
                for key in count:
                    count[key] = count[key] / len(kmers)
            #code = [name, label]
            code = []
            for j in range(2, len(header)):
                if header[j] in count:
                    code.append(count[header[j]])
                else:
                    code.append(0)
            encoding.append(code)
    return encoding

###
def DNC(fastas, **kw):
    base = 'ACGU'

    encodings = []
    dinucleotides = [n1 + n2 for n1 in base for n2 in base]
    header = ['#', 'label'] + dinucleotides

    #encodings.append(header)

    AADict = {}
    for i in range(len(base)):
        AADict[base[i]] = i

    for i in fastas:
        #name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        sequence, label = i[0], i[1]
        #code = [name, label]
        code =[]
        tmpCode = [0] * 16
        for j in range(len(sequence) - 2 + 1):
            tmpCode[AADict[sequence[j]] * 4 + AADict[sequence[j+1]]] = tmpCode[AADict[sequence[j]] * 4 + AADict[sequence[j+1]]] +1
        #print(f'tmpCode {tmpCode}')
        
        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]
        code = code + tmpCode
        #print(f'code {code}')

        encodings.append(code)
    return encodings
    
###
def TNC(fastas, **kw):
    AA = 'ACGU'
    encodings = []
    triPeptides = [aa1 + aa2 + aa3 for aa1 in AA for aa2 in AA for aa3 in AA]
    header = ['#', 'label'] + triPeptides

    #encodings.append(header)

    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i

    for i in fastas:
        #name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        sequence, label = i[0], i[1]
        #code = [name, label]
        code = []
        tmpCode = [0] * 64
        for j in range(len(sequence) - 3 + 1):
            tmpCode[AADict[sequence[j]] * 16 + AADict[sequence[j+1]]*4 + AADict[sequence[j+2]]] = tmpCode[AADict[sequence[j]] * 16 + AADict[sequence[j+1]]*4 + AADict[sequence[j+2]]] +1
        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]
        code = code + tmpCode
        encodings.append(code)
    return encodings
    
###
def ENAC(fastas, window=5, **kw):
    if check_fasta_with_equal_length == False:
        print('Error: for "ENAC" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    if window < 1:
        print('Error: the sliding window should be greater than zero' + '\n\n')
        return 0

    if get_min_sequence_length(fastas) < window:
        print('Error: all the sequence length should be larger than the sliding window :' + str(window) + '\n\n')
        return 0

    AA = kw['order'] if kw['order'] != None else 'ACGU'
    encodings = []
    header = ['#', 'label']

    for w in range(1, len(fastas[0][0]) - window + 2): ###
        for aa in AA:
            header.append('SW.' + str(w) + '.' + aa)
    #encodings.append(header)

    for i in fastas:
        #name, sequence, label = i[0], i[1], i[2]
        sequence, label = i[0], i[1]
        #code = [name, label]
        code = []
        for j in range(len(sequence)):
            if j < len(sequence) and j + window <= len(sequence):
                count = Counter(sequence[j:j + window])
                for key in count:
                    count[key] = count[key] / len(sequence[j:j + window])
                for aa in AA:
                    code.append(count[aa])
        encodings.append(code)
    return encodings
    
###

def binary(fastas, **kw):
    if check_fasta_with_equal_length == False:
        print('Error: for "BINARY" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    AA = 'ACGU'
    encodings = []
    header = ['#', 'label']
    for i in range(1, len(fastas[0][0]) * 4 + 1): ### fastas[0][1] => fastas[0][0]
        header.append('BINARY.F'+str(i))
    #encodings.append(header)

    for i in fastas:
        sequence, label = i[0], i[1]
        #code = [name, label]
        code = []
        for aa in sequence:
            if aa == '-':
                code = code + [0, 0, 0, 0]
                continue
            for aa1 in AA:
                tag = 1 if aa == aa1 else 0
                code.append(tag)
        encodings.append(code)
    return encodings     


###
def CKSNAP(fastas, gap=5, **kw):
    if gap < 0:
        print('Error: the gap should be equal or greater than zero' + '\n\n')
        return 0

    if get_min_sequence_length(fastas) < gap + 2:
        print('Error: all the sequence length should be larger than the (gap value) + 2 = ' + str(gap + 2) + '\n\n')
        return 0

    AA = kw['order'] if kw['order'] != None else 'ACGU'
    encodings = []
    aaPairs = []
    for aa1 in AA:
        for aa2 in AA:
            aaPairs.append(aa1 + aa2)

    header = ['#', 'label']

    for g in range(gap + 1):
        for aa in aaPairs:
            header.append(aa + '.gap' + str(g))
    #encodings.append(header)

    for i in fastas:
        #name, sequence, label = i[0], i[1], i[2]
        sequence, label = i[0], i[1]
        #code = [name, label]
        code = []
        for g in range(gap + 1):
            myDict = {}
            for pair in aaPairs:
                myDict[pair] = 0
            sum = 0
            for index1 in range(len(sequence)):
                index2 = index1 + g + 1
                if index1 < len(sequence) and index2 < len(sequence) and sequence[index1] in AA and sequence[index2] in AA:
                    myDict[sequence[index1] + sequence[index2]] = myDict[sequence[index1] + sequence[index2]] + 1
                    sum = sum + 1
            for pair in aaPairs:
                code.append(myDict[pair] / sum)
        encodings.append(code)
    return encodings
   
###
chemical_property = {
    'A': [1, 1, 1],
    'C': [0, 1, 0],
    'G': [1, 0, 0],
    'T': [0, 0, 1],
    'U': [0, 0, 1],
    '-': [0, 0, 0],
}

def NCP(fastas, **kw):
    if check_fasta_with_equal_length == False:
        print('Error: for "NCP" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    AA = 'ACGU'
    encodings = []
    header = ['#', 'label']

    for i in range(1, len(fastas[0][0]) * 3 + 1): ### fastas[0][1] => fastas[0][0]
        header.append('NCP.F'+str(i))
    #encodings.append(header)

    for i in fastas:
        #name, sequence, label = i[0], i[1], i[2]
        sequence, label = i[0], i[1]
        #code = [name, label]
        code = []
        for aa in sequence:
            code = code + chemical_property.get(aa, [0, 0, 0])
        encodings.append(code)
    return encodings

###
def ANF(fastas, **kw):
    if check_fasta_with_equal_length == False:
        print('Error: for "ANF" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    AA = 'ACGU'
    encodings = []
    header = ['#', 'label']

    for i in range(1, len(fastas[0][0]) + 1): ### fastas[0][1] => fastas[0][0]
        header.append('ANF.' + str(i))
    #encodings.append(header)

    for i in fastas:
        #name, sequence, label = i[0], i[1], i[2]
        sequence, label = i[0], i[1]
        #code = [name, label]
        code = []
        for j in range(len(sequence)):
            code.append(sequence[0: j + 1].count(sequence[j]) / (j + 1))
        encodings.append(code)
    return encodings


###
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

def CalculateMatrix(data, order):
    matrix = np.zeros((len(data[0]) - 2, 64))
    for i in range(len(data[0]) - 2): # position
        for j in range(len(data)):
            if re.search('-', data[j][i:i+3]):
                pass
            else:
                matrix[i][order[data[j][i:i+3]]] += 1
    return matrix


def PSTNPss(fastas, **kw):
    if check_sequences.check_fasta_with_equal_length == False:
        print('Error: for "PSTNP" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    for i in fastas:
        if re.search('[^ACGU-]', i[0]):
            print('Error: illegal character included in the fasta sequences, only the "ACGT[U]" are allowed by this encoding scheme.')
            return 0

    encodings = []
    header = ['#', 'label']
    for pos in range(len(fastas[0][0])-2): ###
        header.append('Pos.%d' %(pos+1))
    #encodings.append(header)

    positive = []
    negative = []
    positive_key = []
    negative_key = []
    for i in fastas:
        #if i[3] == 'training':
            #print(str(i[1]))    
            if str(i[1]) == '1': ###
                #print('OK')
                positive.append(i[0])
                positive_key.append(i[0])
            else:
                negative.append(i[0])
                negative_key.append(i[0])

    nucleotides = ['A', 'C', 'G', 'U']
    trinucleotides = [n1 + n2 + n3 for n1 in nucleotides for n2 in nucleotides for n3 in nucleotides]
    order = {}
    for i in range(len(trinucleotides)):
        order[trinucleotides[i]] = i

    matrix_po = CalculateMatrix(positive, order)
    matrix_ne = CalculateMatrix(negative, order)

    positive_number = len(positive)
    negative_number = len(negative)

    for i in fastas:
        #if i[3] == 'testing':
            #name, sequence, label = i[0], i[1], i[2]
            sequence, label = i[0], i[1]
            #code = [name, label]
            code = []
            for j in range(len(sequence) - 2):
                if re.search('-', sequence[j: j+3]):
                    code.append(0)
                else:
                    p_num, n_num = positive_number, negative_number
                    po_number = matrix_po[j][order[sequence[j: j+3]]]
                    if i[0] in positive_key and po_number > 0:
                        po_number -= 1
                        p_num -= 1
                    ne_number = matrix_ne[j][order[sequence[j: j+3]]]
                    if i[0] in negative_key and ne_number > 0:
                        ne_number -= 1
                        n_num -= 1
                    code.append(po_number/p_num - ne_number/n_num)
                    # print(sequence[j: j+3], order[sequence[j: j+3]], po_number, p_num, ne_number, n_num)
            encodings.append(code)
    return encodings



