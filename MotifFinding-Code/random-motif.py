#!/usr/bin/python
import numpy as np
import random
import itertools as iter
from scipy import stats

def generateRand(L):
    out = ''
    for i in range(L):
        out += np.random.choice(['a','c','g','t'])
    return out

def mutateDeterministic(s,m):
    #m = number of mutations
    if m == 0:
        return s
    places = [i for i in range(len(s))]
    random.shuffle(places)
    locs = places[:m]
    for loc in locs:
        alph = ['a','c','g','t']
        alph.remove(s[loc]) #remove so that the mutation always gives something different
        s = s[:loc] + np.random.choice(alph) + s[loc+1:]
    return s
    
def mutateRandom(s,p):
    #p = probability of mutation at any given site
    for i in range(len(s)):
        if random.random() < p:
            alph = ['a','c','g','t']
            alph.remove(s[i])
            s = s[:i]+np.random.choice(alph)+s[i+1:]
    return s
    
def insert(s,input,i):
    assert i < len(s), 'Oops, insertion location out of bounds.'
    return s[:i]+input+s[i:]

def insertRand(s,input):
    return insert(s,input,np.random.randint(len(s)))

def insertRandMutated(s,input,type='prob',p=0.5,m=1):
    assert type == 'prob' or type == 'det', 'Allowed types are prob and det.'
    if type == 'prob':
        return insertRand(s, mutateRandom(input,p))
    if type == 'det':
        return insertRand(s, mutateDet(input,m))
        

def profile(motifs):
    #motifs is a list of kmers (k=len(motifs[0]))
    #Output is a profile matrix where each column is a distribution
    #on 4 outcomes.
    
    dict = {'a': 0, 'c': 1, 'g': 2, 't': 3}
    
    #using laplace
    profile = np.ones( (4,len(motifs[0])) )
    for row in motifs:
        for i in range(len(row)):
            profile[dict[row[i]],i] += 1
    return profile/(len(motifs)+2)
    
def sampleFromProfile(profile):
    output = ''
    for i in range(len(profile[0])):
        output += random.choices(['a','c','g','t'], profile[:,i])[0] 
        #because the output of choice is a list (with a single element)
    return output

def likeliestString(profile):
    #outputs the most likely string
    dict = {0: 'a', 1: 'c', 2: 'g', 3: 't'}
    
    s = ''
    for i in range(len(profile[0])):
        s+= dict[np.argmax(profile[:,i])]
    return s


# prof = profile(['aaaaa','ccccc','ggggg','ttttt'])
# 
# print( prof )    
# print( prof[:,0] )   
# for i in range(10):
#     print(sampleFromProfile(prof))
# 


def generateMotifs(data,k):
    motifs = [data[i][np.random.randint(len(data[i])-k):] for i in range(len(data))]
    motifs = [row[:k] for row in motifs]
    return motifs

def complement(s1):
    s2 = s1[::-1]
    s2.replace('A','B')
    s2.replace('T','A')
    s2.replace('B','T')
    
    s2.replace('C','B')
    s2.replace('G','C')
    s2.replace('B','G')
    return s2
    
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

def score(motifs,data):
    return sum([Hamming(motifs[i],data[i]) for i in range(len(data))])
    
def profileProb(kmer,prof):
    #returns the probability of a kmer according to profile
    
    assert len(prof[0]) == len(kmer), 'Oops, profile size doesnt match k.'
    
    dict = {'a': 0, 'c': 1, 'g': 2, 't': 3}
    prob = 1
    for i in range(len(kmer)):
        prob *= prof[ dict[kmer[i]] ][i]
    return prob
    
def newMotifsDet(data,prof):
    k = len(prof[0])
    probMatrix = [[profileProb(row[i:i+k],prof) for i in range(len(row)-k)] for row in data]
    maxes = [np.argmax(matRow) for matRow in probMatrix]
    return [data[i][maxes[i]:maxes[i]+k] for i in range(len(maxes))]
    
def newMotifsProb(data,prof):
    k = len(prof[0])
    probMatrix = [[profileProb(row[i:i+k],prof) for i in range(len(row)-k)] for row in data]
    sampledPositions = [random.choices([i for i in range(len(probRow))], probRow) for probRow in probMatrix]
    return [data[i][sampledPositions[i][0]:sampledPositions[i][0]+k] for i in range(len(data))]


def motifFinderDet(data,k,stillIterations):
    motifs = generateMotifs(data,k)
    stopParam = 0
    while True:
        prof = profile(motifs)
        sc = score(motifs,data)
        motifsAlt = newMotifsDet(data,prof)
        scAlt = score(motifsAlt,data)
        if scAlt <= sc:
            motifs = motifsAlt
            if scAlt < sc:
                stopParam = 0
            else:
                stopParam += 1
        else:
            stopParam += 1
        if stopParam == stillIterations:
            break
        print(stopParam,end='\r')
    return motifs
    
def motifFinderProb(data,k,stillIterations):
    motifs = generateMotifs(data,k)
    stopParam = 0
    while True:
        prof = profile(motifs)
        sc = score(motifs,data)
        motifsAlt = newMotifsProb(data,prof)
        scAlt = score(motifsAlt,data)
        if scAlt <= sc:
            motifs = motifsAlt
            if scAlt < sc:
                stopParam = 0
            else:
                stopParam += 1
                # print(stopParam)
        else:
            stopParam += 1
            # print(stopParam)
        if stopParam == stillIterations:
            break
        # print(stopParam,end='\r')
    return motifs

def findConsensus(motifs):
    dict = {'a': 0, 'c': 1, 'g': 2, 't': 3}
    dict2 = {0: 'a', 1: 'c', 2: 'g', 3: 't'}
    consensus = ''
    for i in range(len(motifs[0])):
        n = np.zeros(4)
        for j in range(len(motifs)):
            n[dict[motifs[j][i]]] += 1
        c = dict2[np.argmax(n)]
        consensus += c
    return consensus
    
    


def instantiate(seqLength,seqN,motifLength,pMutate,iter):
    realMotif = motifLength*'a'
    data = [insertRandMutated(generateRand(seqLength),realMotif,'prob',pMutate) for i in range(seqN)]
    motifs = motifFinderProb(data,motifLength,iter)
    prof = profile(motifs)
    entropy = 0
    entropy = sum([stats.entropy([row[i] for row in prof], base=2) for i in range(len(prof[0])) ])    
    #this is bullshit, show me entropy!
    #also, calculating the *bits of entropy
    consensus = findConsensus(motifs)
    
    return [motifs,consensus,entropy]



# motifLength = 20
# pMutate = 0.001
# iter = 2    

# for x in range(1,2):
#     seqN = x*50
#     f = open(f"N{seqN}.txt","w+")
#     f.write(f"#Motif Length = {motifLength}, p mutation = {pMutate}, iterations = {iter}\n")
#     f.write(f"#Number of Seqs = {seqN}\n")
#     f.write("#Seq Length, Entropy\n")
#     for y in range(1,14):
#         avgEntr = 0
#         avgN = 10
#         for z in range(avgN):
#             seqLength= y*50
#             instance = instantiate(seqLength,seqN,motifLength,pMutate,iter)
#             print("\n n of seq = ", seqN,", length = ", seqLength,"\n")
#             # print("\n".join(instance[0]), "\n")
#             print("consensus = ", instance[1])
#             print("entropy = ", instance[2])
#             print("distance = ", Hamming('a'*motifLength, instance[1]))
#             avgEntr += instance[2]/avgN
#             f.write(f"{seqLength}, {instance[2]}\n") #have all the data points in plot, not only avg
#         print(f"Avg Entropy, for Length = {seqLength}, is = {avgEntr}")
#     f.close()
# 
    
f = open("upstream250.txt")
dataCaps = []
while True:
    aux = f.readline()
    if aux == '':
        break
    dataCaps.append(f.readline())

aux = None
f.close()
data=[]
uncap = {'A': 'a', 'C': 'c', 'G': 'g', 'T': 't'}
for datum in dataCaps:
    x = ''
    for a in datum:
        if a=='\n':
            break
        x += uncap[a]
    data.append(x)

k = 18

for i in range(3):
    iterations = input("number of iterations = ")

    motifs = motifFinderProb(data,k,iterations)
    prof = profile(motifs)
    entropy = 0
    entropy = sum([stats.entropy([row[i] for row in prof], base=2) for i in range(len(prof[0])) ])    
    consensus = findConsensus(motifs)
    # print(f"-----------\nk = {k}\n----------\n")
    # print(motifs)
    print(consensus,", ",entropy)
    # print(entropy)




# stillIterations = 2
# datapts = 300
# realMotif = 'aaaaaaaaaaaaaaaa'
# k = len(realMotif)
# data = [insertRandMutated(generateRand(400),realMotif,'prob',p=0.001) for i in range(datapts)]
# print("\n".join(data),"\n")
# motif = motifFinderProb(data,k,stillIterations)
# print("\n".join(motif),"\n")
# print(score(motif,data))


# 
# def motifFinderProb(data,k):
#     motifs = generateMotifs(data,k)
#     stopParam = 0
#     while True:
#         prof = profile(motifs)
#         sc1 = score(motifs,data)
#         motifsAlt = [sampleFromProfile(prof) for l in range(len(data))]
#         if (score(motifsAlt,data) <= sc1):# =<
#             motifs = motifsAlt
#             stopParam = 0
#         else:
#             stopParam += 1
#         if stopParam == 100:
#             break
#     return motifs
# 
# def motifFinderDet(data,k):
#     motifs = generateMotifs(data,k)
#     stopParam = 0
#     while True:
#         prof = profile(motifs)
#         sc1 = score(motifs,data)
#         motifChoice = likeliestString(prof)
#         if score(motifChoice,data) <= sc1:
#             motifs = motifsAlt
#             if score(motifChoice,data) < sc1:
#                 stopParam = 0
#             else:
#                 stopParam += 1
#         else:
#             stopParam += 1
#         if stopParam == 100:
#             break
#     return motifs
# 
# datapts = 40
# realMotif = 'aaaaaaaaaa'
# k = len(realMotif)
# data = [insertRandMutated(generateRand(20),realMotif,'prob',p=0.001) for i in range(datapts)]
# print("\n".join(data),"\n")
# motif = motifFinderProb(data,k)
# print("\n".join(motif),"\n")
# print(profile(motif))
# 





    