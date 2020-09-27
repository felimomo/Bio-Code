#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt

#OriC is the origin of coding (where the genome transcription starts).
#Also, because of the directionality of transcription, there is one side
#that lives longer as a single strand (higher mutation rate) --this is the
#part upstream of the origin in terms of transcription (opposite direction
#than the DNA direction). Since C has a higher mutation rate, there is a 
#C deficit wrt G in this region.
#
#
#On the downstream part, the opposite is true: there is a C surplus wrt G (since
#its complement is a downstream part with C deficit). 
#
#
#Following the direction of the genome (which would be upstream of the coding direction)
#we find that #G-#C increases up to the origin point, and then decreases. (Up to 
#statistical fluctuations.)
#
#I'll start by plotting this difference function

f=open("../Genome.txt")
gen = f.read()
diff=0
plotx = []
ploty = []
for x in gen:
    if x=='G':
        diff +=1
    if x=='C':
        diff -=1
    plotx.append(gen.index(x))
    ploty.append(diff)
    
plt.scatter(ploty,plotx)
plt.show()



