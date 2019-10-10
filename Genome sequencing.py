# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 23:19:19 2018

@author: Shufan Chen
"""
import numpy as np


data = open("E:\\download\\dataset_6207_2.txt","r")
Text = data.readlines()
text = []
for each in Text:
    text.append(each.rstrip())
text1 = Text[0].rstrip()
text2 = Text[1].rstrip()
data.close()



def Composition(Text,k):
    result = []
    for i in range(len(Text)-k+1):
        result.append(Text[i:i+k])
    return result

a = Composition('0100011101',3)
 
def GenomePath(Strings):
    k=len(Strings[0])
    l = len(Strings)+k-1
    result = ""
    for each in Strings[::k]:
        result = result+each
    n = l%k
    if n==0:
        return result
    else:
        result=result+Strings[-1][-n:]
        return result         
GenomePath(['GGC', 'GCT', 'CTT', 'TTA', 'TAC', 'ACC', 'CCA','CAT'])



def Overlap_Graph(Strings):
    temp={}
    result = []
    for each in Strings:
        temp[each[1:]] = [each]
        for i in Strings:
            if each[1:] == i[:-1]:
                temp[each[1:]].append(i)
    for link in temp.values():
        if len(link) == 1:
            continue
        else:
            t = link[0]+" -> "+",".join(link[1:])
            result.append(t)
    return(result)
Overlap_Graph(['ATGCG','GCATG','CATGC','AGGCA','GGCAT','GGCAC'])



def DeBruijn_link(String,k):
    strings = []
    for i in range(len(String)-k+2):
        strings.append(String[i:i+k-1])
    result = {}
    for e in range(len(strings)-1):
        if strings[e] in result.keys():
            result[strings[e]].append(strings[e+1])
        else:
            result[strings[e]] = [strings[e+1]]
    return result

DeBruijn_link('TAATGCCATGGGATGTT',3)
c=[]
for k,v in zip(a.keys(),a.values()):
    t = k+" -> "+",".join(v[:])
    c.append(t)

def DeBruijn_linkk(Strings):
    result = {}
    for each in Strings:
        if each[:-1] in result.keys():
            result[each[:-1]].append(each[1:])
        else:
            result[each[:-1]] = [each[1:]]
    return result

DeBruijn_linkk(['GAGG','CAGG','GGGG','GGGA','CAGG','AGGG','GGAG'])
    
#Eulerian_Circuit_Algorithm
def ECA_overall(text,ppt):
    graph = {}
    if len(ppt) == 2:
        graph[ppt[1]]=[]
    for each in text:
        temp = each.split(" -> ")
        k=temp[0]
        v=temp[1].split(",")
        graph[k]=v
    circuit = []
    i = ppt[0]
    def ECA(i):
        if not graph[i]:
            circuit.insert(0,i)
        else:
            while graph[i]:
                j = graph[i].pop(0)
                ECA(j)
            circuit.insert(0,i)
        return "->".join(circuit)
    return ECA(i)
    
    

#Eulerian_Path_Algorithm
Strings=['0 -> 2','1 -> 3','2 -> 1','3 -> 0,4','6 -> 3,7','7 -> 8','8 -> 9','9 -> 6']
def find_thread(Strings): 
    Out = {}
    for link in Strings:
        t = link.split(" -> ")
        Out[t[0]]=0
    for each in Strings:
        temp = each.split(" -> ")
        v=temp[1].split(",")
        Out[temp[0]] += len(v)
        for e in v:
            if e in Out.keys():
                Out[e] += -1
            else:
                Out[e] = -1
    for k,v in Out.items():
        if v == 1:
            k1 = k
        elif v==-1:
            k2 = k
    return [k1,k2]
ppt = find_thread(text)

ECA_overall(text,ppt)

#DNA reconstruction
def DR(Strings):
    Out = {}
    d = DeBruijn_linkk(Strings)
    for link in d:
        Out[link]=0
    for i,j in d.items():
        Out[i] += len(j)
        for e in j:
            if e in Out.keys():
                Out[e] += -1
            else:
                Out[e] = -1
    ppt=[]
    for k,v in Out.items():
        if v == 1:
            ppt.insert(0,k)
        elif v == -1:
            ppt.insert(1,k)
    if ppt:
        d[ppt[1]]=[]
    else:
        ppt.append(Strings[0][:-1])  
    circuit = []
    i = ppt[0]
    def ECA(i):
        if not d[i]:
            circuit.insert(0,i)
        else:
            while d[i]:
                j = d[i].pop(0)
                ECA(j)
            circuit.insert(0,i)
        return circuit
    return GenomePath(ECA(i))

DR(['AAAT','AATG','ACCC','ACGC','ATAC','ATCA','ATGC','CAAA','CACC','CATA','CATC','CCAG','CCCA','CGCT','CTCA','GCAT','GCTC','TACG','TCAC','TCAT','TGCA'])
DR(['000','001','010','011','100','101','110','111'])

a = DR(text)

fileObject = open('sample1.txt', 'w')  
for n in c:  
    n = str(n)
    fileObject.write(n)  
    fileObject.write('\n')  
fileObject.close()


#de_bruijn sequancing
def de_bruijn(k, n):
    """
    de Bruijn sequence for alphabet k
    and subsequences of length n.
    """
    try:
        # let's see if k can be cast to an integer;
        # if so, make our alphabet a list
        _ = int(k)
        alphabet = list(map(str, range(k)))

    except (ValueError, TypeError):
        alphabet = k
        k = len(k)

    a = [0] * k * n
    sequence = []

    def db(t, p):
        if t > n:
            if n % p == 0:
                sequence.extend(a[1:p + 1])
        else:
            a[t] = a[t - p]
            db(t + 1, p)
            for j in range(a[t - p] + 1, k):
                a[t] = j
                db(t + 1, t)
    db(1, 1)
    return "".join(alphabet[i] for i in sequence)

de_bruijn(2,9)



#Reconstructing a String from the Paired de Bruijn Graph
import networkx as nx

#(3,1)
text = ['TAA|GCC','AAT|CCA','ATG|CAT','TGC|ATG','GCC|TGG','CCA|GGG','CAT|GGA','ATG|GAT','TGG|ATG','GGG|TGT','GGA|GTT']
#(3,2)
text = ['AAT|CAT','ATG|ATG','ATG|ATG','CAT|GAT','CCA|GGA','GCC|GGG','GGG|GTT','TAA|CCA','TGC|TGG','TGG|TGT']
#(4,2)
text = ['GACC|GCGC','ACCG|CGCC','CCGA|GCCG','CGAG|CCGG','GAGC|CGGA']
#(3,1)
text = ['ACC|ATA','ACT|ATT','ATA|TGA','ATT|TGA','CAC|GAT','CCG|TAC','CGA|ACT','CTG|AGC','CTG|TTC','GAA|CTT','GAT|CTG','GAT|CTG','TAC|GAT','TCT|AAG','TGA|GCT','TGA|TCT','TTC|GAA']

def Graph(array, k):
    graph = {}
    l = len(text)
    for string in array:
        prefix = string[0:k-1] + string[k] + string[k+1:2*k]
        suffix = string[1:k] + string[k] + string[k+2:]
        try:
            graph[prefix].append(suffix)
        except(KeyError):
            graph[prefix] = [suffix]    
    return graph,l

def EulerPath(text,k,d):
    g,l= Graph(text, k)
    G = nx.MultiDiGraph(g)
    nodes = list(G.nodes())
    print(len(nodes))
    for node in nodes:
        if G.in_degree(node) < G.out_degree(node):
           path = list(nx.edge_dfs(G, node))
           break
    string = path[0][0][:k-1]
    for i in range (0, len(path)):
        string += path[i][1][k-2]
    for j in range(l-k-d, len(path)):
        string += path[j][1][-1]
    return string

EulerPath(text,3,1)
c = Graph(text,120)

# Maximal Non-Branching Paths in a Graph
def Graph_1(array):
    graph = {}
    for string in array:
        prefix = string[:-1]
        suffix = string[1:] 
        try:
            graph[prefix].append(suffix)
        except(KeyError):
            graph[prefix] = [suffix]    
    return graph
text = ['ATG','ATG','TGT','TGG','CAT','GAT','GGA','AGA']
G = Graph_1(text)
values = []
for val in G.values():
    values.extend(val)
def MaximalNonBranchingPaths(Graph):
    paths = []
    values = []
    for val in Graph.values():
        values.extend(val)
    nodes = [n for n in Graph if len(Graph[n]) != values.count(n) or values.count(n) > 1]
    for node in nodes:
        for edge in Graph[node]:
            path = node+edge[-1]
            while edge in Graph and values.count(edge)==1 and len(Graph[edge])==1:
                edge = Graph[edge][0]
                path = path + edge[-1]
            paths.append(path)
    return paths
            
a = MaximalNonBranchingPaths(G)   
c=[]
for each in a:
    print(each)

text = ['1 -> 2,3,5','2 -> 1,4','3 -> 2,5','4 -> 1,2,5','5 -> 3']  
G={}  
for link in text:
    t = link.split(" -> ")    
    if len(t[1]) > 1:
        G[t[0]]=t[1].split(",")
    else:
        G[t[0]]=[t[1]]
        
def MaximalNonBranchingPaths_r(Graph):
    paths = []
    values = []
    delit = []
    for val in Graph.values():
        values.extend(val)
    nodes = [n for n in Graph if len(Graph[n]) != values.count(n) or values.count(n) > 1]
    for node in nodes:
        for edge in Graph[node]:
            path = node+"->"+edge
            delit.append(edge)
            while edge in Graph and values.count(edge)==1 and len(Graph[edge])==1:
                edge = Graph[edge][0]
                delit.append(edge)
                path = path + "->" + edge
            paths.append(path)
    for rm in delit:
        values.remove(rm)
    for e in values:
        path = e 
        while Graph[e]:
            path = path + "->" + Graph[e][0]
            e = Graph[e].pop()
            values.remove(e)
        paths.append(path)
    return paths        
a = MaximalNonBranchingPaths_r(G)
for each in a:
    print(each)

