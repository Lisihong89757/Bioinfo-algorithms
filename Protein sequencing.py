# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 14:21:17 2018

@author: Shufan Chen
"""

data = open("E:\\download\\dataset_104_8.txt","r")
Text = data.readlines()
text = []
for each in Text:
    text.append(each.rstrip())
text1 = Text[0].rstrip()
text2 = Text[1].rstrip()
data.close()




from itertools import product

ref = [''.join(x) for x in product("ACGU", repeat=3)]
PRO = 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV Y YSSSS CWCLFLF'
RNA_PRO = {}
for i in range(len(ref)):
    RNA_PRO[ref[i]] = PRO[i]

def translation(text):
    l = len(text)
    protein = ""
    for i in range(0,l,3):
        rna = text[i:i+3]
        temp = RNA_PRO[rna]
        protein = protein + temp
    protein = protein.rstrip()
    return protein


translation('CCCCGUACGGAGAUGAAA')
translation(text)


RNA_PRO
#Peptide Encoding Problem
trans = {"A":"U","T":"A","G":"C","C":"G"}
def peptide_encode(text,protein):
    Pro = []
    rna_complement = ""
    for each in text:
        rna_complement = rna_complement+trans[each]
    rna = text.replace("T","U")
    rna_complement = rna_complement[::-1]
    pro = ""
    pro_com = ""
    for i in range(len(rna)-2):
        pro = pro+RNA_PRO[rna[i:i+3]]
        pro_com = pro_com + RNA_PRO[rna_complement[i:i+3]]
    k = len(protein)
    for i in range(len(pro)-(k-1)*3):
        if pro[i:(k-1)*3+i+1:3] == protein:
            Pro.append(text[i:i+k*3])
        if pro_com[i:(k-1)*3+i+1:3] == protein:
            Pro.append(text[-i-k*3:-i])
    return Pro

peptide_encode(text[0],"HWSPYSTFGT")


#Cyclospectrum cyclic peptide


mass = ['G57','A71','S87','P97','V99','T101','C103','I113','L113','N114','D115','K128','Q128','E129','M131','H137','F147','R156','Y163','W186']
ProMass = {}
for each in mass:
    ProMass[each[0]]=int(each[1:])

def extended_mass_table():
    extended_list = {}
    for i in range(57,201):
        extended_list[chr(i)] = i
    return(extended_list)
extended_mass_table()

 
c = 'LEQN'

def Cyclospectrum(text):
    all_pro = []
    for i in range(len(text)):
        n=text[i:]+text[:i]
        all_pro.append(n[:-1])
    result = []
    for each in all_pro:
        def count(seq):
            if len(seq) == 1:
                result.append(ProMass[seq])
                return ProMass[seq]
            else:
                wgh = ProMass[seq[-1]]+count(seq[:-1])
                result.append(wgh)
                return wgh 
        count(each)
    result.sort()
    result.append(sum(result[:len(text)]))
    result.insert(0,0)
    return result

o= [str(i) for i in Cyclospectrum('QCV')]
" ".join(o)




#Cyclospectrum cyclic peptide not recursive
def peptide_masses(peptide):
    '''
    convert peptite string to a list of masses
    '''
    global ProMass
    return list(map(lambda k:ProMass[k],list(peptide)))

c = peptide_masses('VAQ')

def peptide_mass_spectrum(pmass, cyclic = True):
    ''' 
    convert list of peptide masses to spectrum
    '''
    s = [0, ]
    ll = list(pmass)    
    n = len(ll)
    it = None
    if cyclic:
        ll.extend(pmass[:-1])
        s.append(sum(pmass))
        it = [(i,j) for i in range(n) for j in range (i+1,i+n)]
    else:
        it = [(i,j) for i in range(n) for j in range (i+1,n+1)]
        
    for (i,j) in it:
            subpeptide_mass = sum(ll[i:j])
            s.append(subpeptide_mass)
    
    return sorted(s)

peptide_mass_spectrum(c,False)


# BFCyclopeptideSequencing
def counting_peptides_with_given_mass(mass):
    '''
    compute the number of peptides of given total mass.
    '''
    aam = sorted(list(set(ProMass.values())), reverse = True)
    md = {0:1}
    for i in range(min(aam), mass+1):
        for m in aam:
            if i-m in md:
                md[i] = md[i-m] + md.get(i,0)
    return md[mass]

counting_peptides_with_given_mass(1000)




#Cyclopeptide  Sequencing

def spectrum_consistent(p,s, cyclic = False):
    lsp = peptide_mass_spectrum(p, cyclic = cyclic)
    for i in lsp:
        if not i in s:
            return False
    return True

def cyclopeptide_sequencing(spectrum):
    '''
    find all peptides consistent with a given spectrum
    '''
    lp = [[]]
    res =  []
    lmass = list(set(extended_list.values()))
    spectrum.sort(reverse = True)
    parent_mass = max(spectrum)
    def expand(a):
        exp = []
        for i in a:
            for j in lmass:
                p = list(i)
                p.append(j)
                exp.append(p)
        return exp
    while lp:
        lp = expand(lp)
        for p in list(lp):
            if sum(p) == parent_mass:
                if spectrum_consistent(p, spectrum, cyclic = True):
                    res.append(p)
                lp.remove(p)
            elif not spectrum_consistent(p, spectrum):
                lp.remove(p)
    return res


a = [0,71,99,101,103,128,129,199,200,204,227,230,231,298,303,328,330,332,333]

for each in cyclopeptide_sequencing(a):
    print("-".join([str(i) for i in each]))
    
    
peptide_mass_spectrum(peptide_masses("VYYEVDWTMGRQIDPDEYPIAQCTRHRATILTLPDWQM"))    
#Cyclopeptide Scoring Problem    


def scoring(text,spectrum,cyclic=True):
    if (type(text) == str):
        comp = peptide_mass_spectrum(peptide_masses(text),cyclic=cyclic)
    else:
        comp = peptide_mass_spectrum(text,cyclic=cyclic)
    comp.sort()
    spectrum.sort()
    score = 0
    i=0
    j=0
    while i < len(comp) and j < len(spectrum):
        if(spectrum[j]==comp[i]):
            j += 1
            i += 1
            score += 1
        elif(spectrum[j]>comp[i]):
            i+=1
        else:
            j+=1
    return score
    
       
spec = [0,99,113,114,128,227,257,299,355,356,370,371,484]    
scoring("VVCLKRQLIRPFDFWIMHQIAIIPHASDPFDDEHQFPDELH",text, False) 


#Trim
def leaderboard_trim(leaderboard, spectrum, N):
    '''
    output the N highest-scoring linear peptides on Leaderboard lp
    with respect to Spectrum
    '''
    # need for trimming ?
    if len(leaderboard)<=N:
        return leaderboard
    
    #build a dict of peptide:score
    d = {tuple(e): scoring(e, spectrum, cyclic = False) for e in leaderboard }
    #排序
    ll = sorted(d, key=d.get, reverse = True)
    min_score = d[ll[N-1]]
    tp = []
    for e in ll:
        if (d[e]>=min_score):
            tp.append(list(e))
        else:
            # cut off loop optimization
            return tp
    return tp
leaderboard_trim(a,b,6)


#Spectral Convolution Problem
def spectral_convolution(spectrum):
    spectrum.sort()
    conv = []
    n = len(spectrum)
    for i in range(n):
        for j in range(i+1,n):
            diff = spectrum[j]-spectrum[i]
            if diff > 0:
                conv.append(diff)
    return conv

import operator
from collections import Counter

def leaderboard_cyclopeptide_sequencing(spectrum, N, M = 20, convolution = False):
    '''
    find all peptides approximatively consistent 
    with a given spectrum
    '''
    lp = [[]]
    top = []
    top_score = 0
    spectrum.sort(reverse = True)
    lmass = None
    if not convolution:
        lmass = list(set(ProMass.values()))
    else:
        #build a specific alphabet from the convoluted spectrum
        candidate = [k for k in spectral_convolution(spectrum) if k>= 57 and k <= 200]
        dconv = Counter(candidate)
        lconv = sorted(dconv.items(), key=operator.itemgetter(1), reverse = True)
        min_freq = 0
        lmass = []
        for (k, v) in lconv:
            if (M > 0):
                lmass.append(k)
                min_freq = v
                M -= 1
            elif (M == 0) and (v == (min_freq)):
                #handle ties
                lmass.append(k)
    parent_mass = max(spectrum)
    
    def expand(a):
        exp = []
        for i in a:
                for j in lmass:
                    p = list(i)
                    p.append(j)
                    exp.append(p)
        return exp
    
    while lp:
        lp = expand(lp)
        for p in list(lp):
            if sum(p) == parent_mass:
                p_score = scoring(p, spectrum)
                if p_score >= top_score:
                    top_score = p_score
                    top = p
                    print(p,p_score)
            elif sum(p) > parent_mass:
                lp.remove(p)
        lp = leaderboard_trim(lp, spectrum, N)
    return top

a = [0,71,113,129,147,200,218,260,313,331,347,389,460]
text = text[1].split()
text = [int(i) for i in text]
a = leaderboard_cyclopeptide_sequencing(spectrum, N=1000,M=20,convolution=True)
a = [str(i) for i in a]
print("-".join(a))

scoring([128, 99, 163, 128, 114, 147, 186, 194, 163], spectrum)

scoring('MAMA',[0,71,98,99,131,202,202,202,202,202,299,333,333,333,503])
scoring("PEEP",[0,97,129,129,129,194,226,323,323,355,452],cyclic=False)
spetral = []
spectrum = [0,57,118,179,236,240,301]
for i, v in enumerate(spectrum):
    for s in spectrum[:i]:
        if v != s:
            spetral.append(v-s)
spetral.sort()
filtered = [c for c in spetral if c >= 57 and c <= 200]
        
elements = sorted([(c, filtered.count(c)) for c in set(filtered)],
                   key=lambda tup: tup[1], reverse=True)
