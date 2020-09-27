#!usr/bin/python
import networkx as nx
import numpy as np
import sys
sys.path.insert(0,'../')
import Tools.hamming as ham

def different(kmer1,kmer2,pMut):
    #identifies the two kmers as coming from mutations of different kmers.
    #For this I want to use some tail bounds and set a low error prob,
    #but for now just do simple stuff
    assert len(kmer1)==len(kmer2)
    k = len(kmer1)
    dist = ham.hamming(kmer1,kmer2)
    # print((k+np.sqrt(k))*pMut)
    if dist <= (k+np.sqrt(k))*pMut:
        return False
    else:
        return True
        
def edgeBool(node1,node2,kmers,pMut):
    #should 2 nodes be connected by a directed edge?
    for kmer in kmers:
        if (not different(node1,kmer[:-1],pMut)) and (not different(node2,kmer[1:],pMut)):
            return True
    # if !(different(node1,kmers[1:])) and !(different(node2,kmers[:-1])):
    #     return True
    #not this one to preserve direction!
    else:
        return False

def deBrujn(kmers,pMut):
    deBrujn = nx.DiGraph()
    for kmer in kmers:
        pref = kmer[:-1]
        suf = kmer[1:]
        # print(pref,"\n",suf,"\n")
        if all([different(pref,node,pMut) for node in list(deBrujn.nodes())]):
            deBrujn.add_node(pref)
        if all([different(suf,node,pMut) for node in list(deBrujn.nodes())]):
            deBrujn.add_node(suf)
    print([len(vert) for vert in deBrujn.nodes()])
    edges = [(node1,node2) for node1 in list(deBrujn.nodes()) for node2 in list(deBrujn.nodes()) if edgeBool(node1,node2,kmers,pMut)]
    # print(edges)
    deBrujn.add_edges_from(edges)
    return deBrujn
        
    
    
