#!/usr/bin/python

def reverse(s):
    return "".join([s[-1-i] for i in range(len(s))])

def complement(s):
    exch = {'a':'t', 't':'a', 'g':'c', 'c':'g'}
    sRev = reverse(s)
    return "".join([exch[c] for c in sRev])

def hamming(s1,s2):
    assert len(s1)==len(s2)
    l = len(s1)
    h1 = sum([1 for i in range(l) if s1[i]!=s2[i]])
    h2 = sum([1 for i in range(l) if s1[i]!=complement(s2)[i]])
    return min(h1,h2)