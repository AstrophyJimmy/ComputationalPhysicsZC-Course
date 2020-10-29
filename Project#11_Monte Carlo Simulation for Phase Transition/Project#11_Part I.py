"""
Assingment No. 11 Part I
Name: Mohamed Gamal Zaid
ID: 201700399
"""
import numpy as np 
from numpy import exp as E
import matplotlib.pyplot as plt 
from matplotlib.pyplot import cm



J=1
T=4
H=0
n=10
total = np.power(n,2)
ts=1100
nCut = 100
plot = False


def interactingSpinsIndices(i,j,s=1):
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


def normalize(m,n):
    return m/np.power(n,2)

def MontCarlo(T, H, n, ts):
    spins = np.ones([n,n])
    Ms = np.zeros(ts) #Magnetization
    Es = np.zeros(ts) #Average Energies
    deltaEs = np.zeros(ts) #Delta Energy = <E^2> - <E>^2 = AvgEsSq - AvgSqEs
    for t in range(ts):
        #First I remove the boundary spins to allow for looping without worrying about the BCs
        summ = 0
        summSq = 0
        for i in range(n):
            for j in range(n):
                inds = interactingSpinsIndices(i,j)
                if (t!=0):
                    Eflip = 2*J*np.sum( [np.product([spins[i,j], spins[ind1,ind2]]) for ind1,ind2 in inds[1:]])
                    if Eflip <= 0:
                        spins[i,j]=-spins[i,j]
                    elif Eflip > 0:
                        r=np.random.rand()
                        BoltzFactor = E(-Eflip/T)
    
                        if(r <= BoltzFactor):
                            spins[i,j]=-spins[i,j]
                #calculating energies per each iteration per each spin
                #E = -Eflip/4
                AvgEng = -0.5*J*np.sum( [np.product([spins[i,j], spins[ind1,ind2]]) for ind1,ind2 in inds[1:]])
                summ += AvgEng
                summSq += np.power(AvgEng,2)
        if plot:
            plt.matshow(spins)
            plt.savefig("run\\"+str(t)+".jpeg")
            plt.close("all")
        #if(not Ms[0]): print(spins)
        Ms[t]=np.sum(spins)/total

        Es[t]=summ/total
        deltaEs[t] = summSq/total - np.power(Es[t],2)
        #AvgEsSq[t]=summSq/total
        #AvgSqEs[t]=np.power(Ms[t],2)
    #return Ms, Es, AvgEsSq, AvgSqEs

    return Ms, Es, deltaEs
        
#M,Energy,dEnergy=MontCarlo(T, H, n, ts)

Temp = np.arange(0.5,5,0.05)
#Temp = [2.25]
nPoints = np.size(Temp)
color=cm.rainbow(np.linspace(0, 1, nPoints))
Mag = np.zeros(nPoints)
AvgEnergy = np.zeros(nPoints)
specificHeat = np.zeros(nPoints)
import time
st = time.perf_counter()
for i,T in enumerate(Temp):
    st1 = time.perf_counter()
    T=np.round(T,2)
    print("Run #"+str(i+1)+" of "+str(nPoints))
    #M, Energy, AvgE_Sq, AvgSq_E = MontCarlo(T, H, n, ts)
    M, Energy, dEnergySq= MontCarlo(T, H, n, ts)
    
    Mag[i] = np.mean(M[nCut:])
    AvgEnergy[i] = np.mean(Energy[nCut:])
    #specificHeat[i] = (np.mean(AvgE_Sq[nCut:]) - np.mean(AvgSq_E[nCut:]))/np.power(T,2)
    specificHeat[i] = np.mean(dEnergySq[nCut:])*10/(np.power(T,2))
    Title = "Figure 8.6"
    if T in [1.0, 1.5, 2.25, 4.0]:
        plt.plot(range(ts),M,c=color[i],label='T='+str(T),alpha=0.8)
    en1 = time.perf_counter()
    print("Time of run #"+str(i+1)+" took: "+str(np.round((-st1+en1)/60,3))+" mins")

print("Total Time is: "+ str(np.round((-st+en1)/60,3))+" mins")

plt.ylabel("Magnetization")
plt.xlabel("Time")
plt.title(Title)
plt.grid(alpha=0.2)  
plt.legend()
plt.savefig(Title+".jpeg",dpi=300,pad_inches=0.5)
plt.close("all") 


plt.scatter(Temp, Mag, alpha=0.7, label="10x10 Lattice")
plt.xlabel("Temprature")
plt.ylabel("Magnetization")
Title = "Figure 8.7"
plt.title(Title)   
plt.legend() 
plt.grid(alpha=0.2)  
plt.savefig(Title+".jpeg",dpi=300,pad_inches=0.5)
plt.close("all") 


plt.scatter(Temp, AvgEnergy, alpha=0.7, label="10x10 Lattice")
plt.xlabel("Temprature")
plt.ylabel("Average Energy")
Title = "Figure 8.8"
plt.title(Title)   
plt.legend() 
plt.grid(alpha=0.2)  
plt.savefig(Title+".jpeg",dpi=300,pad_inches=0.5)
plt.close("all") 


plt.scatter(Temp, specificHeat, alpha=0.7, label="10x10 Lattice")
plt.xlabel("Temprature")
plt.ylabel("Specific Heat")
Title = "Figure 8.9"
plt.title(Title)   
plt.legend() 
plt.grid(alpha=0.2)  
#plt.savefig(Title+".jpeg",dpi=300,bbox_inches='tight',pad_inches=1)
plt.savefig(Title+".jpeg",dpi=300,pad_inches=0.5)
plt.close("all") 