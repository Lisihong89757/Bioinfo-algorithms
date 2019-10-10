# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 21:34:56 2018

@author: Shufan Chen
"""
def hammingdis(text1,text2):
    distance = 0
    for i in range(len(text1)):
        if text1[i] != text2[i]:
            distance += 1
    return distance

def Neighbors(pattern, d):
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {"A","C","G","T"}
    neighborhood = set([])
    pre_nei = Neighbors(pattern[1:], d)
    for each in pre_nei:
        if hammingdis(pattern[1:], each) < d:
            for ele in ["A","C","G","T"]:
                neighborhood.add(ele + each)
        else:
            neighborhood.add(pattern[0] + each)
    return neighborhood

#brute force search
def ModifEnumeration(Dna, k, d):
    temp = []
    for i in range(len(Dna)):
        temp.append([])
    for each_s in Dna:
        n=Dna.index(each_s)
        for i_k in range(len(each_s)-k+1):
            p = Neighbors(each_s[i_k:i_k+k],d)
            temp[n].extend(list(p))
    Patterns = set(temp[0]).intersection(*temp[1:])
    return list(Patterns)

Dna = ['ATTTGGC',
'TGCCTTA',
'CGGTATC',
'GAAAATT']
a = ModifEnumeration(Dna,3,1)


# compute entropy
import numpy as np
import pandas as pd
import math
b = ["TCGGGGGTTTTT",
    "CCGGTGACTTAC",
    "ACGGGGATTTTC",
    "TTGGGGACTTTT",
    "AAGGGGACTTCC",
    "TTGGGGACTTCC",
    "TCGGGGATTCAT",
    "TCGGGGATTCCT",
    "TAGGGGAACTAC",
    "TCGGGTATAACC"]
c = [[],[],[],[],[],[],[],[],[],[]]
for i in range(len(b)):
    for each in b[i]:
        c[i].append(each)
c = np.array(c)
c = pd.DataFrame(c)


df_new = pd.DataFrame(columns = ["A","C","G","T"])
for i in range(12):
    temp = c[i].value_counts()
    df_new.loc[i] = temp/10
df_new = df_new.fillna(0)

result = 0
for each in df_new:
    for x in df_new[each]:
        if x != 0:
            t = x*math.log2(x)
            result += t
#Meddianstring
def DBPAS(Pattern, Dna):
    k = len(Pattern)
    result_d = 0
    for each in Dna:
        HD = 10000
        for i in range(len(each)-k+1):
            d = hammingdis(each[i:i+k], Pattern)
            if HD > d:
                HD = d
        result_d += HD
    return result_d
DBPAS("AAA",['TTACCTTAAC', 'GATATCTGTC', 'ACGGCGTTCG', 'CCCTAAAGAG', 'CGTCAGAGGT'])

from itertools import product
def MedianString(Dna,k):
    dtn = 100000
    ref = [''.join(x) for x in product("ACGT", repeat=k)]
    for each in ref:
        d = DBPAS(each, Dna)
        if dtn >= d:
            dtn = d
            Median = each
    return Median
MedianString(['CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC',
'GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC',
'GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG'],7)


data = open("E:\\download\\dataset_160_9.txt","r")
Text = data.readlines()
Pattern = Text[0].rstrip()
Array = []
for each in Text[1:]:
    h = each.rstrip()
    Array.append(h)
Array = np.array(Array)

Dna = Text[1].split(" ")
data.close()
DBPAS(Pattern, Dna)

k = int(Pattern)
Dna= Text[1:11]
MedianString(Dna,k)

#Greedy Modif search
def PPkP(Text, k, Array):
    s = ["A","C","G","T"]
    comp = 0
    result = Text[0:k]
    for idx in range(len(Text)-k):
        num = 1
        f = Text[idx:idx+k]
        for each,j in zip(f, range(k)):
            i = s.index(each)
            p = Array[i,j]
            num = num*p
        if num > comp:
            comp = num
            result = f
    return result

PPkP("ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT",5)

def scoring(Modifs):
    k = len(Modifs[0])
    n = len(Modifs) 
    s = ["A","C","G","T"]
    score = 0
    array = np.zeros(shape=(4,k))
    
    for m in Modifs:
        for each,j in zip(m,range(k)):
            i = s.index(each)
            array[i,j] += 1
    for i in range(k):
        score = score + n - max(array[:,i])
    return score
    
# greeding modif search with laplace rule of succession   
def GreedyModifSearch(Dna,k,t):
    Best_m = []
    s = ["A","C","G","T"]
    for each_s in Dna:
        Best_m.append(each_s[0:k])
    for i in range(len(Dna[0])-k+1):
        Modif = []
        arr = np.zeros(shape=(4,k))
        fst = Dna[0][i:i+k]
        for j in range(len(fst)):
            arr[s.index(fst[j]),j] += 1
        Modif.append(fst)
        n=1
        for each in Dna[1:]:
            #laplace rule of succession
            arr_t = (arr+1)/(n+4)
            temp = PPkP(each,k,arr_t)
            Modif.append(temp)
            for i in range(len(temp)):
                arr[s.index(temp[i]),i] += 1
            n += 1
        if scoring(Best_m) > scoring(Modif):
            Best_m = Modif
    return Best_m

Dna = ['GGCGTTCAGGCA','AAGAATCAGTCA','CAAGGAGTTCGC','CACGTCAATCAC','CAATAATATTCG']
            
a = GreedyModifSearch(Array,12,25)       
                    
PPkP('AAGAATCAGTCA', 3, Array)        
Array = np.array([[0,0,0],
                  [0,0,1],
                  [1,1,0],
                  [0,0,0]])    





fileObject = open('sample1.txt', 'w')  
for n in a:  
    n = str(n)
    fileObject.write(n)  
    fileObject.write('\n')  
fileObject.close()            
            