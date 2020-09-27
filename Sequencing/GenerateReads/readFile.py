#!/usr/bin/python

def simpleReader(fileName):
    #reads file with only complete genome, returns genome as a lower-case string
    f = open("../InFiles/"+fileName)
    S = f.read().replace('\n','')
    if S[0].isupper():
        SLower = ''
        dict = {'A': 'a', 'C': 'c', 'G': 'g', 'T': 't'}
        for s in S:
            SLower += dict[s]
        return SLower
    return S
    

    