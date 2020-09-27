#!usr/bin/python
import GenerateReads.generateGenome as gen
import MakeDeBrujn.deBrujn as graph
import GenerateReads.sampleReads as samp
import Reconstruct.reconstruct as rec
import networkx as nx
import matplotlib.pyplot as plt
genome = gen.literGen('Shakespeare1.txt',70)
kmers = samp.kmerSamples(genome, 25, 0.01, 250)
deB = graph.deBrujn(kmers,0.01)
print(nx.is_connected(deB.to_undirected()),"\n")
print(nx.adjacency_matrix(deB))
eu = nx.eulerian_path(deB)
euler = list(eu)
# print(euler)
# print(euler[0])
# print(euler[0][0])
# print(euler[0][0][0])
genRec = rec.naiveReconstruct(euler)
print(genome,"\n")
print(genRec)