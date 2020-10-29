"""
Assingment No. 11 Part II
Name: Mohamed Gamal Zaid
ID: 201700399
"""
import numpy as np 
from numpy import exp as E
import matplotlib.pyplot as plt 
from matplotlib.pyplot import cm
import time


J=1
T=5
H=0
n=20
total = np.power(n,2)
ts=1100
nCut = 100
plot = False
plotCorr = False
s=1


def interactingSpinsIndices(i,j,s):
    """
    Find all possible combination between the indices of the 
    desired spins and all of its neighbors
    to calculate Eflip Accurately
    and allowing for Boundary Conditions
    indices ordered as follows: 0:center, 1:left, 2:right, 3:up, 4:down
    each time 0:i, 1:j
    """
    indices=np.array([(i,j),(i,j-s),(i,j+s),(i-s,j),(i+s,j)],dtype=int)
    
    #We have four corners and four edges at corners we need two indices
    #at edges we need just one
    right = n-j-s-1
    down = n-i-s-1
    left = j-s
    up = i-s
    
    if left<0: #left edge 1
        indices[1,1] = (n+left) #j
    
    elif right<0: #right edge 2
        indices[2,1] = (-right-1) #j
    
    if up<0: #upper edge 3
        indices[3,0] = (n+up) #i
    
    elif down<0: #lower edge 4
        indices[4,0] = (-down-1) #i
    return indices


def MontCarloCorrelation(T, H, n, ts, s):
    spins = np.ones([n,n])
    correlFunc = np.zeros(ts)
    for t in range(ts):
        #First I remove the boundary spins to allow for looping without worrying about the BCs
        for i in range(n):
            for j in range(n):
                inds = interactingSpinsIndices(i,j,1)
                indsSpaced = interactingSpinsIndices(i,j,s)
                
                Eflip = 2*J*np.sum( [np.product([spins[i,j], spins[ind1,ind2]]) for ind1,ind2 in inds[1:]])
                if t>=nCut:
                    correlFunc[t] += np.mean( [np.product([spins[i,j], spins[ind1,ind2]]) for ind1,ind2 in indsSpaced[1:]])

                if Eflip <= 0:
                    spins[i,j]=-spins[i,j]
                elif Eflip > 0:
                    r=np.random.rand()
                    BoltzFactor = E(-Eflip/T)
                    if(r <= BoltzFactor):
                        spins[i,j]=-spins[i,j]

        if plot:
            plt.matshow(spins)
            plt.savefig("run\\"+str(t)+".jpeg")
            plt.close("all")

        if plotCorr:
            plt.matshow(spins)
            plt.savefig("runCorr\\"+str(s)+"_"+str(t)+".jpeg")
            plt.close("all")        

        
        correlFunc=correlFunc/total
    return correlFunc
        
st = time.perf_counter()

nSteps = np.int(n/2)
steps = np.linspace(1,nSteps,nSteps)
TempRange = [1.5, 2, 2.25, 3, 5]
color=cm.rainbow(np.linspace(0, 1, len(TempRange)))


for q, T in enumerate(TempRange):    
    
    f_i=np.zeros(nSteps)
    st1 = time.perf_counter()

    for i,s in enumerate(steps):    
        print("for T= "+ str(T)+" & s= " + str(s))

        f_i[i]=np.sum(MontCarloCorrelation(T, H, n, ts, s))    
    plt.plot(steps ,f_i ,c=color[q],label='T='+str(T),alpha=0.8)
    plt.scatter(steps ,f_i ,c=color[q].reshape(1,4),alpha=0.8)
    en1 = time.perf_counter()
    print("It took: "+str(np.round((en1-st1)/60,3))+" Mins")
print("In total it took: "+str(np.round((en1-st)/60,3))+" Mins")

plt.xlabel("i=Distance [Steps]")
plt.ylabel(r"$f(i)=<s_i s_0>$",fontsize = 12)
Title = "Figure 8.10"
plt.title(Title)   
plt.legend() 
plt.grid(alpha=0.2)
plt.ylim(-0.1,1)  
plt.savefig(Title+".jpeg",dpi=300,pad_inches=0.5)
#plt.close("all") 