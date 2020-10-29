"""
Assingment No. 11 Part III
Name: Mohamed Gamal Zaid
ID: 201700399
"""
import numpy as np 
from numpy import exp as E
import matplotlib.pyplot as plt 
from matplotlib.pyplot import cm
import time



J=1
T=1
H=-5
n=10
total = np.power(n,2)
ts=1100
nCut = 100
plot = False



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
    #print(indices)
    return indices


def MontCarlo(T, H, n, ts, sign=1):
    spins = sign *np.ones([n,n])
    Ms = np.zeros(ts) #Magnetization

    for t in range(ts):
        #First I remove the boundary spins to allow for looping without worrying about the BCs
        
        for i in range(n):
            for j in range(n):
                inds = interactingSpinsIndices(i,j,s=1)
                if (t!=0):
                    Eflip = 2*(J*np.sum( [np.product([spins[i,j], spins[ind1,ind2]]) for ind1,ind2 in inds[1:]])+spins[i,j]*H)
                    if Eflip <= 0:
                        spins[i,j]=-spins[i,j]
                    elif Eflip > 0:
                        r=np.random.rand()
                        BoltzFactor = E(-Eflip/T)
    
                        if(r <= BoltzFactor):
                            spins[i,j]=-spins[i,j]

        if plot:
            plt.matshow(spins,cmap = cm.viridis)
            plt.savefig("run\\"+str(t)+".jpeg")
            plt.close("all")
        Ms[t]=np.sum(spins)/total

    return Ms


st = time.perf_counter()

TempRange = [1, 2, 2.5, 3, 5]
color=cm.rainbow(np.linspace(0, 1, len(TempRange)))
Hs = np.arange(-5,5.5,0.5)
nH = np.size(Hs)


for q, T in enumerate(TempRange):    
    
    M=np.zeros(nH)
    st1 = time.perf_counter()

    for i,H in enumerate(Hs):    
        print("for T= "+ str(T)+" & H= " + str(H))
        M[i]=np.mean(MontCarlo(T, H, n, ts))    
        
    plt.plot(Hs ,M ,c=color[q],label='T='+str(T),alpha=0.7)
    plt.scatter(Hs ,M ,c=color[q].reshape(1,4),alpha=0.7)
    en1 = time.perf_counter()
    print("It took: "+str(np.round((en1-st1)/60,3))+" Mins")

print("In total it took: "+str(np.round((en1-st)/60,3))+" Mins")


Title = "Figure 8.11, 8.12"
plt.ylabel("M")
plt.xlabel("H")
plt.title(Title)
plt.grid(alpha=0.2)  
plt.legend()
plt.savefig(Title+".jpeg",dpi=300,pad_inches=0.5)
#plt.close("all") 

