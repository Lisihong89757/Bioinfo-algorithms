# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 13:29:59 2018

@author: Shufan Chen
"""
from itertools import product

def ComputingFrequencies(Text,k):
    ind_num = {}
    ref = [''.join(x) for x in product("ACGT", repeat=k)]
    for pattern in ref:
        ind_num[pattern] = 0
    for i in range(len(Text)-(k-1)):
        Pattern = Text[i:i+k]
        ind_num[Pattern] += 1
    return ind_num


def Clumpfinder(Genome,k,L,t):
    Ltclump = set([])
    ref_dict = ComputingFrequencies(Genome[0:L],k)
    for i,j in ref_dict.items():
        if j == t:
            Ltclump.add(i)
    for i in range(len(Genome)-L):
        ref_dict[Genome[i:i+k]] -= 1
        ref_dict[Genome[i+L-k:i+L]] += 1
        if ref_dict[Genome[i+L-k:i+L]] == t:
            Ltclump.add(Genome[i+L-k:i+L])
    return list(Ltclump)

data = open("E:\\download\\sample.txt","r")
Text = data.readlines()
Text1 = Text[0].rstrip()
k = int(Text[1].rstrip())
data.close()

fileObject = open('sample1.txt', 'w')  
for n in a:  
    n = str(n)
    fileObject.write(n)  
    fileObject.write('\n')  
fileObject.close()





dict1 = {"A":0,"C":1,"G":2,"T":3}
def P2N(Pattern):
    if len(Pattern) == 0:
        return 0
    value = dict1[Pattern[-1]]
    return 4*P2N(Pattern[:-1])+value


dict2 = {0:"A",1:"C",2:"G",3:"T"}
result = ""
def N2P(Index, k):
    if k == 1:
        dna = dict2[Index]
        return dna
    temp = Index%4
    Index = Index//4
    dna = dict2[temp]
    return N2P(Index,k-1)+ dna




