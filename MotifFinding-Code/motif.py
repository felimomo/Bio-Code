#!/usr/bin/python
import numpy as np
import itertools as iter
import os

def complement(s1):
    s2 = s1[::-1]
    s2.replace('A','B')
    s2.replace('T','A')
    s2.replace('B','T')
    
    s2.replace('C','B')
    s2.replace('G','C')
    s2.replace('B','G')
    return s2

def iseq(s1,s2):
    #take into account complement possibility
    if s1==s2 or complement(s1)==s2:
        return True
        
    return False

def Hamming(s1,s2):
    if len(s1) > len(s2):
        s3 = s1
        s1 = s2
        s2 = s3
        s3 = None
    #always have it so that len(s1) =< len(s2)
    
    ham = len(s1)
    hamcomp = len(s1)
    if len(s1)==len(s2):
        for i in range(len(s1)):
            if s1[i] == s2[i]:
                ham -= 1
            if s1[i] == complement(s2)[i]:
                hamcomp -= 1
        return min([ham,hamcomp])
    
    for i in range(len(s2)-len(s1)):
        if Hamming(s1, s2[i:i+len(s1)]) < ham:
            ham = Hamming(s1, s2[i:i+len(s1)])
    return ham

#use for finding musical motifs accross different songs later!!
#but in music motifs can happen more than once in a song

def ConsensusFind(data,k,alph):
    #data will be a list of len(data)-many strings (possibly with different lengths)
    #alph is the alphabet
    distance = k*len(data)
    consensus = []
    iteration = 0
    for kmer in iter.product(alph,repeat=k):
        iteration += 1
        distHypoth = sum([Hamming(kmer,datum) for datum in data]) # this could be expensive to compute so only do it once
        if distHypoth < distance:
            distance = distHypoth
            consensus = kmer
        print(f"consensus={consensus}, it.={iteration}/{4**k}, dist={distHypoth} .. {distance}", end='\r') 
    return [consensus, distance]

def ConsensusToMotif(consensus,data):
    assert len(consensus)<min(map(len, data)), "Oh no! Data too short, consensus too long!"
    motifs = []
    for datum in data:
        position = np.argmin([ Hamming(consensus, datum[i:i+len(consensus)]) for i in range(len(datum)-len(consensus)) ])
        motifs.append(datum[position:position+len(consensus)])    
    return motifs
    
def MotifFind(data,k,alph):
    consensus = ConsensusFind(data,k,alph)
    return [consensus,ConsensusToMotif(consensus[0],data)]
    
f = open('../Genome.txt')
length = 100
n = 10
alph = ['A', 'C', 'G', 'T']
data = []
k = 6
for i in range(n):
    start =np.random.randint(0, os.path.getsize('../Genome.txt')-length)
    f.seek(start)
    data.append(f.read(length))

motifs = MotifFind(data,k,alph)
print("\n",motifs)
            
# data = [ [1,2,3,4,5,6,7,8,9,10,0,0,0,0,1,0,1], [10,7,1,0,6,0,0,4,0,0,8,4,3,8,2,0,4],[7,8,4,7,5,6,2,1,0,0,0,0,0,7,9,5,7]]
# k=6
# alph = [i for i in range(10)]
# 
# reachingCon = ConsensusFind(data,k,alph)
# print(reachingCon)
# motif = MotifFind(reachingCon[0],data)
# print(motif)

            
            
            
            
            
            
            
    