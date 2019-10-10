# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 13:35:18 2018

@author: Shufan Chen
"""
import numpy as np
import pandas as pd

data = open("E://download//2BreakDistance.txt","r")
Text = data.readlines()
text = []
for each in Text:
    text.append(each.rstrip())
text1 = Text[0].rstrip()
text2 = Text[1].rstrip()
data.close()

def DPChange(money, coins):
    """
    Using dynamic programming method to calculate the minimum coins for 
    change
    """
    Mincoins = [0]
    for m in range(1, money+1):
        Mincoins.append(min([Mincoins[m-i] for i in coins if m >= i]) + 1)
    return Mincoins[money]

DPChange(17564,[19,12,11,9,8,5,3,1])
DPChange(25,[3,2,1])

def ManhattanTourist(n, m, down, right):
    s = np.zeros((n+1,m+1),dtype='int8')
    for i in range(n):
        s[i+1,0] = s[i,0] + down[i,0]
    for j in range(m):
        s[0,j+1] = s[0,j] + right[0,j]
    for i in range(1,n+1):
        for j in range(1,m+1):
            s[i,j] = max(s[i-1,j]+down[i-1,j], s[i,j-1] + right[i,j-1])
    return s[n,m]


down = np.array([[1,0,2,4,3],
                 [4,6,5,2,1],
                 [4,4,5,2,1],
                 [5,6,8,5,3]])
right = np.array([[3,2,4,0],
                 [3,2,4,2],
                 [0,7,3,3],
                 [3,3,0,2],
                 [1,3,2,2]])

n = text.index("-")
down = np.array([[int(i) for i in each.split()] for each in text[1:n]])
right = np.array([[int(i) for i in each.split()] for each in text[n+1:]])
ManhattanTourist(15,13,down,right)



def LCSBackTrack(v,w):
    n,m = (len(v),len(w))
    Backtrack = [[None]*m for _ in range(n)]
    s = np.zeros((n+1,m+1),dtype='int64')
    def sub_cost(i, j):
        if v[i-1] == w[j-1]:
            return 1
        return 0
    for i in range(1,n+1):
        for j in range(1,m+1):
            s[i,j] = max(s[i-1,j], s[i,j-1],s[i-1,j-1]+sub_cost(i,j))
            if s[i,j] == s[i-1,j]:
                Backtrack[i-1][j-1] = "d"
            elif s[i,j] == s[i,j-1]:
                Backtrack[i-1][j-1] = "r"
            elif s[i,j] == s[i-1,j-1] + 1 and v[i-1] == w[j-1]:
                Backtrack[i-1][j-1] = "di"
    return s,Backtrack


c = LCSBackTrack(text[0],text[1])

def outputLCS(backtrack, v):
    LCS = []
    i = len(backtrack)-1
    j = len(backtrack[0])-1
    while i >= 0 and j >= 0:
        if backtrack[i][j] == "di":
            LCS.append(v[i])
            i = i-1
            j = j-1
        elif backtrack[i][j] == "d":
            LCS.append(v[i])
            i = i-1
        else:
            LCS.append("-")
            j = j-1
    LCS = LCS[::-1]
    LCS = "".join(LCS)
    return LCS
                 
outputLCS(c[1],text[0])



#not good enough
def topologypath(sp, ep, graph):
    points = [sp]
    distance = {}
    distance[sp] = 0
    path = [[sp]]
    paths = []
    while points:
        s = []
        temp_path = []
        i = 0
        for point in points:
            if point in graph:
                edges = graph[point]
                for v in edges:
                    n = edges[v] + distance[point]
                    g = path[i].copy()
                    g.append(v)
                    temp_path.append(g)
                    if v not in distance:
                        distance[v] = n
                        s.append(v)
                    else:
                        s.append(v)
                        if distance[v] < n:
                            distance[v] = n
            i += 1
        print(s)
        path = temp_path[:]
        print(path)
        for each in path:
            if ep in each:
                paths.append(each) 
        points = s[:]
    return paths,distance[ep]
        
    
   





graph = {0:{1:7,2:4},2:{3:2},1:{4:1},3:{4:3,5:2}}   
topologypath(8,26,graph)    

def creatgraph(text):
    graph = {}
    for each in text:
        each = each.split("->")
        i = int(each[0])
        edge = [int(n) for n in each[1].split(":")]
        if i not in graph:
            edge = {edge[0]:edge[1]}
            graph[i] = edge.copy()
        else:
            graph[i][edge[0]]=edge[1]
    return graph
            
graph = creatgraph(text)            
  
    
    


#Global Alignment Problem

AAl = text.pop(0)
AAl = AAl.replace(" ","")
text = [each.split() for each in text]
data = np.array(text)
df = np.delete(data, np.s_[0], axis=1)

def AAtrans(seq):
    return [AAl.index(i) if i in AAl else i for i in seq]
 

def LCSBackTrack(v,w,dataf,sigma):
    n,m = (len(v),len(w))
    v = AAtrans(v)
    w = AAtrans(w)
    Backtrack = [[None]*m for _ in range(n)]
    s = np.zeros((n+1,m+1),dtype='int64')
    s[0,:] = [-sigma*i for i in range(m+1)]
    s[:,0] = [-sigma*i for i in range(n+1)]
    for i in range(1,n+1):
        for j in range(1,m+1):
            t, g, h = s[i-1, j]-sigma, s[i, j-1]-sigma, s[i-1,j-1]+int(dataf[v[i-1],w[j-1]])
            s[i,j] = max(t,g,h)
            #s[i,j] = max(0,t,g,h)
            #if s[i,j] == 0:
                #Backtrack[i-1][j-1] = "sk"
            if s[i,j] == t:
                Backtrack[i-1][j-1] = "d"
            elif s[i,j] == g:
                Backtrack[i-1][j-1] = "r"
            elif s[i,j] == h:
                Backtrack[i-1][j-1] = "di"
            
    return s,Backtrack


c = LCSBackTrack(text[0],text[1],df)

c = LCSBackTrack("PLEASANTLY","MEANLY",df,5)


def outputglobal(backtrack, v, w):
    # Quick lambda function to insert indels.
    insert_indel = lambda word, i: word[:i] + '-' + word[i:]

    # Initialize the aligned strings as the input strings.
    v_aligned, w_aligned = v, w

    # Get the position of the highest scoring cell in the matrix and the high score.
    i, j = len(v), len(w)

    # Backtrack to the edge of the matrix starting at the highest scoring cell.
    while i*j != 0:
        if backtrack[i-1][j-1] == "d":
            i -= 1
            w_aligned = insert_indel(w_aligned, j)
        elif backtrack[i-1][j-1] == "r":
            j -= 1
            v_aligned = insert_indel(v_aligned, i)
        else:
            i -= 1
            j -= 1

    # Prepend the necessary preceeding indels to get to (0,0).
    for repeat in range(i):
        w_aligned = insert_indel(w_aligned, 0)
    for repeat in range(j):
        v_aligned = insert_indel(v_aligned, 0)

    return [v_aligned,w_aligned]
    
                 
outputglobal(c[1],text[0],text[1])    

outputglobal(c[1],"PLEASANTLY","MEANLY")


def outputlocal(backtrack, v, w):
    g1 = []
    g2 = []
    m = np.where(c[0] ==c[0].max())
    i = m[0][0]-1
    j = m[1][0]-1
    while i >= 0 and j >= 0:
        if backtrack[i][j] == "di":
            g1.append(v[i])
            g2.append(w[j])
            i = i-1
            j = j-1
        elif backtrack[i][j] == "d":
            g1.append(v[i])
            g2.append("-")
            i = i-1
        elif backtrack[i][j] == "r":
            g1.append("-")
            g2.append(w[j])
            j = j-1
        else:
            break
    g1 = g1[::-1]
    g2 = g2[::-1]
    g1 = "".join(g1)
    g2 = "".join(g2)
    return [g1,g2]    
    
outputlocal(c[1],text[0],text[1])    
    
    
#Solve the Alignment with Affine Gap Penalties Problem    
    
def LCSBackTrackAG(v,w,dataf,σ,ε):
    n,m = (len(v),len(w))
    Backtrack = [[None]*m for _ in range(n)]
    BU = [[None]*m for _ in range(n)]
    BL = [[None]*m for _ in range(n)]
    md = np.zeros((n+1,m+1),dtype='int64')
    md[1:,0] = [-σ+(i-1)*-ε for i in range(1,n+1)]
    md[0,1:] = [-σ+(i-1)*-ε for i in range(1,m+1)]
    u = np.zeros((n+1,m+1),dtype='int64')
    u[0,1:] = [-σ+(i-1)*-ε for i in range(1,m+1)]
    u[1:,0] = [-σ*2+(i-2)*-ε for i in range(1,n+1)]
    l = np.zeros((n+1,m+1),dtype='int64')
    l[0,1:] = [-σ*2+(i-2)*-ε for i in range(1,m+1)]
    l[1:,0] = [-σ+(i-1)*-ε for i in range(1,n+1)]
    for i in range(1,n+1):
        for j in range(1,m+1):
            t = l[i,j] = max(l[i-1,j]-ε,md[i-1,j]-σ)
            if t == l[i-1,j]-ε:
                BL[i-1][j-1] = "d"
            else:
                BL[i-1][j-1] = "dm"
            
            g = u[i,j] = max(u[i,j-1]-ε,md[i,j-1]-σ)
            if g == u[i,j-1]-ε:
                BU[i-1][j-1] = "r"
            else:
                BU[i-1][j-1] = "rm"
            
            h = md[i-1,j-1]+int(dataf.loc[v[i-1],w[j-1]])
            md[i,j]= max(t,g,h)
            if md[i,j] == h:
                Backtrack[i-1][j-1] = "di" 
            elif md[i,j] == g:
                Backtrack[i-1][j-1] = "r"
            elif md[i,j] == t:
                Backtrack[i-1][j-1] = "d"
             
            
    return md.max(),Backtrack,BU,BL    
c = LCSBackTrackAG(text[0],text[1],df, 11, 1)    


def outputAG(backtrack, bu, bl, v, w):
    g1 = []
    g2 = []
    i = len(backtrack)-1
    j = len(backtrack[0])-1
    while i >= 0 and j >= 0:
        if backtrack[i][j] == "di":
            g1.append(v[i])
            g2.append(w[j])
            i = i-1
            j = j-1
        elif backtrack[i][j] == "d":
            while i >= 0 and j >= 0:
                g1.append(v[i])
                g2.append("-")
                i = i-1
                if bl[i+1][j] == "d":
                    continue
                else:
                    break
        else:
            while i >= 0 and j >= 0:
                g1.append("-")
                g2.append(w[j])
                j = j-1
                if bu[i][j+1] == "r":
                    continue
                else:
                    break
    g1 = g1[::-1]
    g2 = g2[::-1]
    g1 = "".join(g1)
    g2 = "".join(g2)
    return g1,g2    

outputAG(c[0],c[1],c[2],text[0],text[1])



# Middle Edge in Linear Space Problem

def middle_column_score(v, w, scoring_matrix, sigma):
    '''Returns the score of the middle column for the alignment of v and w.'''
    # Initialize the score columns.
    v1 = AAtrans(v)
    w1 = AAtrans(w)
    S = [[i*j*sigma for j in range(-1, 1)] for i in range(len(v)+1)]
    backtrack = [0]*(len(v)+1)
    # Fill in the Score and Backtrack matrices.
    for j in range(1, int(len(w)/2)+1):
        for i in range(0, len(v)+1):
            if i == 0:
                S[i][1] = -j*sigma
            else:
                scores = [S[i-1][0] + int(scoring_matrix[v1[i-1],w1[j-1]]), S[i][0] - sigma, S[i-1][1] - sigma]
                S[i][1] = max(scores)
                backtrack[i] = scores.index(S[i][1])
                
        if j != int(len(w)/2):
            S = [[row[1]]*2 for row in S]
            
    return [row[1] for row in S], backtrack

middle_column_score('YLTNASAELP','YLNSAEM',df,5)

middle_column_score("FPPF","FFPF",df,5)

print()
def middle_edge(v, w, scoring_matrix, sigma):
    '''Returns the middle edge in the alignment graph of v and w.'''

    # Get the score of the middle column from the source to the middle.  The backtrack matrix is unnecessary here.
    source_to_middle = middle_column_score(v, w, scoring_matrix, sigma)[0]
    
    # Get the score of the middle column from the middle to sink.  Reverse the order as the computations are done in the opposite orientation.
    middle_to_sink, backtrack = map(lambda l: l[::-1], middle_column_score(v[::-1], w[::-1]+['', '$'][len(w) % 2 == 1 and len(w) > 1], scoring_matrix, sigma))
    
    # Get the componentwise sum of the middle column scores.
    scores = list(map(sum, zip(source_to_middle, middle_to_sink)))
    
    # Get the position of the maximum score and the next node.
    max_middle = max(range(len(scores)), key=lambda i: scores[i])
    
    if max_middle == len(scores) - 1:
        next_node = (max_middle, int(len(w)/2)+1)
    else:
        next_node = [(max_middle + 1, int(len(w)/2)+1), (max_middle, int(len(w)/2)+1), (max_middle + 1, int(len(w)/2))][backtrack[max_middle]]
    
    return (max_middle, int(len(w)/2)), next_node

middle_edge('TW','WE',df,5)

def space_efficient_global_alignment(v, w, scoring_matrix, sigma):
    '''Return the global alignment of v and w using a linear space algorithm.'''
    def linear_space_alignment2(top, bottom, left, right):
        '''Constructs the global alignment path using linear space.'''
        if left == right:
           return [v[top:bottom], '-'*(bottom - top)]
        elif top == bottom:
           return ['-'*(right - left), w[left:right]]
        elif bottom - top == 1 or right - left == 1:
           c = LCSBackTrack(v[top:bottom], w[left:right], df, 5)
           return outputglobal(c[1],v[top:bottom],w[left:right])
        else:
           # Get the middle edge and the corresponding nodes.
           mid_node, next_node = middle_edge(v[top:bottom], w[left:right], df, 5)
            # Shift the nodes appropriately, as they currently don't alighn with the top/left starting points.
           mid_node = tuple(map(sum, zip(mid_node, [top, left])))
           next_node = tuple(map(sum, zip(next_node, [top, left])))
           
           # Get the character in each alignment corresponding to the current middle edge.
           # (Take the index modulo the string length to avoid IndexErrors if we reach the end of a string but still have -'s to append.)
           current = [['-', v[mid_node[0] % len(v)]][next_node[0] - mid_node[0]], ['-', w[mid_node[1] % len(w)]][next_node[1] - mid_node[1]]]
           
           # Recursively divide and conquer to generate the alignment.
           A = linear_space_alignment2(top, mid_node[0], left, mid_node[1])
           B = linear_space_alignment2(next_node[0], bottom, next_node[1], right)
           
           return [A[i] + current[i] + B[i] for i in range(2)]
    # Get the alignment and alignment score.
    v_aligned, w_aligned = linear_space_alignment2(0, len(v), 0, len(w))
    v = AAtrans(v_aligned)
    w = AAtrans(w_aligned)
    score = sum([-sigma if '-' in pair else int(scoring_matrix[pair]) for pair in zip(v, w)])
    return str(score), v_aligned, w_aligned

space_efficient_global_alignment('PLEASANTLY','MEANLY',df,5)



import copy
#Greedysorting reversal
def greedysort(P):
    n = len(P)
    P = np.array(P)
    result = []
    for i in range(n):
        if np.abs(P[i]) == i+1:
            P[i] = np.abs(P[i])
            result.append(copy.deepcopy(P))
            continue
        if i + 1 in P:
            end = np.where(P == i+1)[0][0]+1
            P[i:end]=P[i:end][::-1]*-1
            result.append(copy.deepcopy(P))
            P[i] = P[i]*-1
            result.append(copy.deepcopy(P))
        else:
            end = np.where(P == -(i+1))[0][0]+1
            P[i:end]=P[i:end][::-1]*-1
            result.append(copy.deepcopy(P))
    return result
    
r = greedysort(w)


def findbreakpoint(P):
    result = 0
    P.append(len(P)+1)
    P.insert(0,0)
    for i in range(len(P)-1):
        if P[i] + 1 != P[i+1]:
            result += 1
    return result

findbreakpoint(w)



def main():
    '''Main call. Reads, runs, and saves problem specific data.'''
    # Read the input data.
    with open('data/textbook/rosalind_5k.txt') as input_data:
        word1, word2 = [line.strip() for line in input_data.readlines()]

    # Get the middle edge.
    middle = middle_edge(word1, word2, BLOSUM62(), 5)

    # Print and save the answer.
    print (' '.join(map(str, middle)))
    with open('output/textbook/Textbook_05K.txt', 'w') as output_data:
        output_data.write(' '.join(map(str, middle)))

if __name__ == '__main__':
    main()

list(df.loc[pair] for pair in zip('PLEASANTLY','MEASNLY'))

sum([-5 if '-' in pair else int(df.loc[pair]) for pair in zip('PL-EASA-NTLY','ME-AS--N--Y')])

def ChromosomeToCycle(Chromosome):
    l = len(Chromosome)
    Nodes = [0]*l*2;
    for i in range(l):
        Chrome = Chromosome[i]
        if Chrome < 0:
            Nodes[2*i] = -Chrome*2;
            Nodes[2*i+1] = -Chrome*2-1;
        else:
            Nodes[2*i] = Chrome*2-1;
            Nodes[2*i+1] = Chrome*2;
    return Nodes
ChromosomeToCycle([4,5,-6])

L = L.split(" ")
L = [int(i) for i in L]
result = ChromosomeToCycle(L)
result = [str(i) for i in result]

def CycleToChromosome(Nodes):
    l = len(Nodes)
    Chromosome = [0]*int(l/2)
    for i in range(int(l/2)):
        if Nodes[2*i] < Nodes[2*i+1]:
            Chromosome[i] = int(Nodes[2*i+1]/2);
        else:
            Chromosome[i] = int(-Nodes[2*i]/2);
    return Chromosome
CycleToChromosome([1,2,4,3,6,5,7,8])   

L = L.split(" ");
L = [int(i) for i in L];
result = CycleToChromosome(L);
result = [str(i) if i < 0 else '+'+str(i) for i in result];
" ".join(result);

def coloredEdge(P):
    Edges = [];
    for chrome in P:
        cycle = ChromosomeToCycle(chrome);
        cycle += [cycle.pop(0)];
        for i in range(0,len(cycle),2):
            Edges.append(cycle[i:i+2]);
    return Edges

L = L[1:-1]
L = L.split(")(")
L_1 = []
for each in L:
    each = each.split(" ")
    L_1.append([int(i) for i in each])
result = coloredEdge(L_1)
result = [tuple(i) for i in result]
result = tuple(result)

def GraphToGenome(Genome):
    temp = [];
    G = [];
    start = Genome[0][0];
    print(start)
    for i in range(len(Genome)-1):
        if abs(Genome[i][1] - start) == 1:
            temp += Genome[i];
            temp = [temp.pop(-1)] + temp;
            G.append(CycleToChromosome(temp));
            temp = [];
            start = Genome[i+1][0];
        else:
            temp += Genome[i];
    temp += Genome[-1];
    temp = [temp.pop(-1)] + temp;
    G.append(CycleToChromosome(temp));
    return G;
GraphToGenome(L)
L = [(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)]

#Calculation of 2-break distance
def twoBreakDistance(P,Q):
    A = coloredEdge(P)
    B = coloredEdge(Q)
    edge_1 = {}
    edge_2 = {}
    for each in A:
        edge_1[each[0]] = each[1];
        edge_1[each[1]] = each[0];
    for each in B:
        edge_2[each[0]] = each[1];
        edge_2[each[1]] = each[0];
    Blocks = len(edge_1)/2;
    def findcycle(a,b,p,lstart,mark = 1):
        if mark == 1:
           if a[p] == lstart:
               del a[p];
               del b[p];
               return 1;
           else:
               p_0 = a[p];
               del a[p];
               del b[p];
               return findcycle(a,b,p = p_0,lstart = lstart,mark = 2)
        else:
            if b[p] == lstart:
                del a[p];
                del b[p];
                return 1;
            else:
                p_0 = b[p];
                del a[p];
                del b[p];
                return findcycle(a,b,p = p_0,lstart = lstart,mark = 1)
    Cycles = 0;
    while len(edge_1):
        Cycles += 1;
        i = list(edge_1.keys())[0];
        findcycle(edge_1,edge_2,i,i);
    
    Distance = Blocks-Cycles;
    return Blocks

P = []             
for each in Text:
    each = each[1:-2]
    each = each.split(")(")
    if len(each) > 1:
        temp = []
        for i in each:
            i = i.split(" ");
            temp.append([int(o) for o in i])
        P.append(temp)
    else:
        each = each.split(" ")
        P.append([int(i) for i in each])
twoBreakDistance(P[0],P[1])


def twoBreakOnGenomeGraph(GenomeGraph,i1,i2,i3,i4):
    for each in GenomeGraph:
        if i1 in each and i2 in each:
            GenomeGraph.remove(each)
            GenomeGraph.append([i1,i3])
        elif i3 in each and i4 in each:
            GenomeGraph.remove(each)
            GenomeGraph.append([i2,i4])
    












