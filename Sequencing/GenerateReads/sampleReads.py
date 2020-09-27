#!/usr/bin/python
import random

def Mutator(s,pMut):
    sOut = ''
    for i in range(len(s)):
        r = random.uniform(0,1)
        if r < pMut:
            alph = 'acgt'
            sOut += random.choice(alph.replace(s[i],''))
        else:
            sOut += s[i]
    return sOut
            

def kmerSample(genome, k, pMut):
    #returns random kmer, possibly mutated (prob pMut)
    length = len(genome)
    strt = random.randint(0,length-k)
    return Mutator(genome[strt:strt+k],pMut)

def kmerSamples(genome, k, pMut, number):
    #returns list of number-many kmers
    return [kmerSample(genome,k,pMut) for i in range(number)]
    
    
