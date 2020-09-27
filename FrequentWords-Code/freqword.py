#!/usr/bin/python
import numpy as np
import itertools

def HighestFreqInWindow(genData,k):
    #genData is the genomic data within the window to consider
    alph = ['A', 'C', 'G', 'T']
    kmers = [[''.join(a)] for a in itertools.product(alph,repeat=k)]
    for x in kmers:
        x.append(0)
    for i in range(len(genData)-k+1):
        for x in kmers:
            if x[0] == genData[i:i+k]:
                x[1] += 1
    return sorted(kmers, key=lambda v: v[1])[-1]
    #returns the most frequent kmer, together with its frequency

k = 6
window = 200
max_len= 3000

f = open("../Genome.txt","r") 
size=len(f.read(max_len))
f.close()
f = open("../Genome.txt","r") 

WindowRecords=[]
position = 0
while position<max_len:
    f.seek(position)
    genData=f.read(window)
    aux = HighestFreqInWindow(genData,k)
    WindowRecords.append([position,aux[0],aux[1]])
    print(f"{position}/"+str(size), end="\r")
    position += 1

f.close()

print( sorted(WindowRecords, key=lambda v: v[2])[-10:-1] )
#sorted outputs the sorted string but doesnt(!!) sort the string and replace the original.

# position = 0
# while True:
#     f.seek(position)
#     word = f.read(k)
#     if len(word)<k:
#         break
#     for x in kmers:
#         if x[0] == word:
#             x[2] += 1
#     position += 1
    
    


# gen = f.read()
# f.close()
# 
# for i in range(len(gen)-k):
#     word = gen[i:i+k]
#     for x in kmers:
#         if x[0] == word:
#             x[2]+=1
# for i in range(25):
#     print(kmers[i])


    
