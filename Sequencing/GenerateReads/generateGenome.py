#!/usr/bin/python
import random

def randGen(length=1000):
    #generate completely random genome sequence
    alph = ['A','B','C','D']
    gen = ''
    for i in range(length):
        gen += random.choice(alph)
    return gen
    
def cleaner(fileName,alphabet,length):
    #returns a string which is fileName but only
    #keeping characters from alphabet. Starting point is random.
    Loc = input('File Location: ')
    f = open(Loc+'/'+fileName)
    xIn = f.read()
    strt = random.randint(0,len(xIn)-length)
    xOut=''
    for x in xIn[strt:strt+length]:
        if any(a==x for a in alphabet):
            xOut += x
    return xOut

def literGen(fileName,length=1000):
    #generate sequence from a text file, whenever a letter is in xList, it adds an x.
    aList = 'abcdef' 
    cList = 'ghijkl' 
    gList = 'mnopqr' 
    tList = 'stuvxy'
    nucleotides = ['a','c','g','t']
    translation = {}
    for n in nucleotides:
        translation = {**translation, **{c:n for c in eval(n+'List')}}
    alph = aList + cList + gList + tList
    cleanText = cleaner(fileName,alph,length)
    
    return "".join([translation[x] for x in cleanText])
    
    

    
    
    
    
    
    
