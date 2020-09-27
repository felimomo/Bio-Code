#!usr/bin/python

def naiveWeave(s1,s2):
    return s1 + s2[-1]
    

def naiveReconstruct(euler):
    #naive reconstruction of genome from eulerian path with no noise
    genome = euler[0][0]
    for i in range(len(euler)):
        genome += euler[i][1][-1]
    return genome