# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 21:54:08 2018

@author: Shufan Chen
"""
Text = "CATTCCAGTACTTCGATGATGGCGTGAAGA"
def minskew(Text):
    num = 0
    result = [0]
    for each in Text:
        if each == "G":
            num += 1
            result.append(num)
        elif each == "C":
            num -= 1
            result.append(num)
        else:
            result.append(0)
    m = min(result)
    Inf=1000000000
    ind=[]
    while min(result) == m:
        ind.append(result.index(m))
        result[result.index(m)] = Inf
    return ind
minskew(Text)        



data = open("E:\\download\\dataset_9_4.txt","r")
Text = data.readlines()
text1 = Text[0].rstrip()
text2 = Text[1].rstrip()
data.close()


text1 = "CTACAGCAATACGATCATATGCGGATCCGCAGTGGCCGGTAGACACACGT"
text2 = "CTACCCCGCTGCTCAATGACCGGGACTAAAGAGGCGAAGATTATGGTGTG"
def hammingdis(text1,text2):
    distance = 0
    for i in range(len(text1)):
        if text1[i] != text2[i]:
            distance += 1
    return distance
hammingdis(text1,text2)

def appro_match(pattern, text, k):
    n = len(pattern)
    result = []
    for i in range(len(text)-n+1):
        if hammingdis(pattern,text[i:i+n]) <= k:
            result.append(i)
    return len(result)

pattern = "CCC"
text= "CATGCCATTCGCATTGTCCCAGTGA"
k=2


appro_match(pattern, text, k)

# similar pattern with mutation less than d
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
a = list(Neighbors("ACGT",3))


#similar pattern with mutation equal to d
def Neighbors_e(pattern, d):
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {"A","C","G","T"}
    neighborhood = set([])
    li = ["A","C","G","T"]
    li.remove(pattern[0])
    pre_nei = Neighbors(pattern[1:], d)
    for each in pre_nei:
        if hammingdis(pattern[1:], each) < d and len(pattern) > len(each)+1:
            for ele in li :
                neighborhood.add(ele + each)
        elif hammingdis(pattern[1:], each) < d and len(pattern) == len(each)+1:
            for ele in li:
                neighborhood.add(ele + each)
        else:
            neighborhood.add(pattern[0] + each)   
    return neighborhood
b = list(Neighbors_e("ACTAGCATAT",2))

#Computing Frequencies With Mismatches
from itertools import product

def ComputingFrequenciesWithMismatches(Text,k,d):
    ind_num = {}
    ref = [''.join(x) for x in product("ACGT", repeat=k)]
    for pattern in ref:
        ind_num[pattern] = 0
    for i in range(len(Text)-(k-1)):
        Pattern = Text[i:i+k]
        neighbor = Neighbors(Pattern, d)
        for each in neighbor:
            ind_num[each] += 1
    return ind_num
ComputingFrequenciesWithMismatches("ACGTTGCATGTCGCATGATGCATGAGAGCT",4,1)


#Computing Frequencies With Mismatches_sorted

def ComputingFrequenciesWithMismatches_s(Text,k,d):
    ind_num = {}
    patterns = []
    for i in range(len(Text)-(k-1)):
        Pattern = Text[i:i+k]
        patterns.extend(list(Neighbors(Pattern, d)))
    patterns_r = set(patterns)
    for each in patterns_r:
        n = patterns.count(each)
        ind_num[each] = n
    m = max(ind_num.values())
    return [i for i,j in ind_num.items() if j == m]
ComputingFrequenciesWithMismatches_s("ACGTTGCATGTCGCATGATGCATGAGAGCT",4,1)
#Computing Frequencies With Mismatches_sorted complement

def ComputingFrequenciesWithMismatches_sc(Text,k,d):
    ind_num = {}
    pair_num = {}
    patterns = []
    for i in range(len(Text)-(k-1)):
        Pattern = Text[i:i+k]
        patterns.extend(list(Neighbors(Pattern, d)))
    patterns_r = set(patterns)
    for each in patterns_r:
        n = patterns.count(each)
        ind_num[each] = n
    li = list(patterns_r)
    for each in li:
        each_c = converse(each)
        if each_c in li:
            li.remove(each_c)
            pair = each + each_c
            num = ind_num[each] + ind_num[each_c]
            pair_num[pair] = num
    m = max(pair_num.values())
    result = [i for i,j in pair_num.items() if j == m]
    t_result = []
    for each in result:
        t_result.append(each[:k])
        t_result.append(each[k:])
    return t_result
ComputingFrequenciesWithMismatches_sc("ACGTTGCATGTCGCATGATGCATGAGAGCT",4,1)



fileObject = open('sample1.txt', 'w')  
for n in a:  
    n = str(n)
    fileObject.write(n)  
    fileObject.write('\n')  
fileObject.close()


def converse(Text):
    Text = Text[::-1]
    Text1 = ""
    for each in Text:
        if each == "A":
            Text1 = Text1+"T"
            continue
        elif each == "T":
            Text1 = Text1+"A"
            continue
        elif each == "G":
            Text1 = Text1+"C"
            continue
        else:
            Text1 = Text1+"G"
    return Text1