# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 22:04:53 2022

@author: Anthony
"""
import scipy as sci
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt


#define variables in terms of a b c



def matrixgen (omega1,omega2,decay1,decay2,gain1,gain2,kappa):
    eigenvals_1 = [] 
    eigenvals_2 = []
    eigenvec_1 = []
    eigenvec_2 = []
    omega0=(omega1+omega2)/2
    delta=omega0-omega1
    a=((-1*delta)/2)+decay1*1j+gain1*1j
    c=((delta)/2)-decay2*1j+gain2*1j
    for i in range(len(a)):
        for k in range(len(c)):
            
                
            abcmatrix= np.array([[a[i]+0j,kappa+0j],[kappa+0j,c[k]+0j]])
    
            eigenvals, eigenvecs = (la.eig(abcmatrix))
                    
           
            eigenvals_1.append(eigenvals[0])
            eigenvals_2.append(eigenvals[1])
            eigenvec_1.append(eigenvecs[0]) 
            eigenvec_2.append(eigenvecs[1])
    return eigenvals_1 ,  eigenvals_2 , eigenvec_1 , eigenvec_2     

    #abcmatrix=[[a[i]+0j,b+0j],[b+0j,c[k]+0j]]
    #return abcmatrix
    
        







def eigenvalue (a , b , c):
    
    evalue1 = 0.5(-1*np.sqrt(a**2-2*a*c+4*a*(b**2)+c**2+0j)+a+c)
    
    evalue2 = 0.5(np.sqrt(a**2-2*a*c+4*a*(b**2)+c**2+0j)+a+c)
        
    return evalue1 , evalue2 

def eigenvector (a , b , c):
    
    x11 = (-1*((-a+c+np.sqrt(a**2+4*(b**2)-2*a*c+c**2)/(2*b))))
    x21 = (-1*((-a+c+np.sqrt(a**2+4*(b**2)-2*a*c+c**2)/(2*b))))
    
    vec1=[ x11 , 1 ]
    vec2=[ x21 , 1 ]
    
    return vec1 , vec2 




x, y, u, v = matrixgen(1, 3, 0.2, 0.3, np.linspace(0,1,9), np.linspace(0,1,9), 4)
print(u , v)

origin=[0,0]
 

plt.quiver(origin, u[1], color=['r'], scale=21)
plt.quiver(origin, v[1], color=['b'], scale=21)
plt.show()

#general eigenvalue/eigenvector solver for a 2x2 a,b,c matrix with wolfram alpha

#method 2: numpy eigenvalue/eigenvector solver
