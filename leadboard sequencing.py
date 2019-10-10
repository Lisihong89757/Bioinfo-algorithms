# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 11:33:32 2018

@author: Shufan Chen
"""
mass = ['G57','A71','S87','P97','V99','T101','C103','I113','L113','N114','D115','K128','Q128','E129','M131','H137','F147','R156','Y163','W186']
ProMass = {}
for each in mass:
    ProMass[each[0]]=int(each[1:])
   
masses = [57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186]

spectrum = [0,71,113,129,147,200,218,260,313,331,347,389,460]    
n=1000

from collections import defaultdict

def cyclopeptide(peptide):
    def subpeptide(peptide, pos, length):
        if pos+length <= len(peptide):
            return peptide[pos: pos+length]
        else:
            return peptide[pos:] + peptide[:length + pos - len(peptide)]
    return [subpeptide(peptide, p, l) for p in range(len(peptide))
            for l in range(1, len(peptide) + 1)]


def cyclospectrum(peptide):
    return sorted([sum(p) for p in cyclopeptide(peptide)] + [0])


def expand_list(peptides, masses):
    if len(peptides) == 0:
        return [[m] for m in masses]
    return [p + [m] for p in peptides for m in masses]


def score(peptide, spectrum_map):
    peptide_map = defaultdict(int)
    for p in peptide:
        peptide_map[p] += 1
    return sum([min(peptide_map[k], spectrum_map[k])
               for k in peptide_map.keys()])


def cut(peptides, spectrum_map, n):
    scores = sorted([(p, score(cyclospectrum(p), spectrum_map))
                    for p in peptides],
                    key=lambda tup: tup[1], reverse=True)
    leaders = [s for s in scores if n >= len(scores) or s[1] >= scores[n][1]]
    return [leader[0] for leader in leaders]


def Convolution_cyclopeptide_sequencing(spectrum,M,N,convolution = False):   
    spectrum_map = defaultdict(int)
    for i in spectrum:
        spectrum_map[i] += 1
    masses = None
    
    #convolution or not
    if not convolution:
        masses = [57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186]
    else:
        spetral = []
        for i, v in enumerate(spectrum):
            for s in spectrum[:i]:
                if v != s:
                    spetral.append(v-s)
        spetral.sort()
        filtered = [c for c in spetral if c >= 57 and c <= 200]
        
        elements = sorted([(c, filtered.count(c)) for c in set(filtered)],
                          key=lambda tup: tup[1], reverse=True)
        
        leaders = [e for e in elements if e[1] >= elements[M][1]]
        
        masses = sorted([l[0] for l in leaders])
    leaders = []
    best = (0, 0)
    for i in range(500):
        expanded = expand_list(leaders, masses)
        leaders = cut([p for p in expanded], spectrum_map, N)
        top_score = score(cyclospectrum(leaders[0]), spectrum_map)
        if top_score >= best[1]:
            best = (leaders[0], top_score)
            print('best(%s): %s' % (best[1], '-'.join(str(i) for i in best[0])))
        leaders = [p for p in leaders if sum(p) <= max(spectrum)]
        if len(leaders) == 0:
            break
    return 'best(%s): %s' % (best[1], '-'.join(str(i) for i in best[0]))
    
    
spectrum = [57,57,71,99,129,137,170,186,194,208,228,265,285,299,307,323,356,364,394,422,493]
Convolution_cyclopeptide_sequencing(spectrum,20,1000, convolution = True)


