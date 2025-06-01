###Calculating and displaying structure functions of neutrino DIS, Fe is given as an example###

from parton import mkPDF
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import integrate as int

#making Iron structure functions for neutrinos#
def FFe(x,Q2):
    dFe=pdfFe.xfxQ2(1,x,Q2)
    dbFe=pdfFe.xfxQ2(-1,x,Q2)
    
    uFe=pdfFe.xfxQ2(2,x,Q2)
    ubFe=pdfFe.xfxQ2(-2,x,Q2)
    
    sFe=pdfFe.xfxQ2(3,x,Q2)
    sbFe=pdfFe.xfxQ2(-3,x,Q2)
    
    cFe=pdfFe.xfxQ2(4,x,Q2)
    cbFe=pdfFe.xfxQ2(-4,x,Q2)
    
    bFe=pdfFe.xfxQ2(5,x,Q2)
    bbFe=pdfFe.xfxQ2(-5,x,Q2)
    
    F2Fe = 2*(dFe + sFe  + ubFe + bFe + cbFe)
    F2bFe= 2*(uFe + cFe + dbFe + sbFe + bbFe)

    F3Fe = (2*(dFe + sFe - ubFe + bFe - cbFe))/x
    F3bFe= (2*(uFe + cFe - dbFe - sbFe - bbFe))/x 
    
    return [F2Fe,F3Fe,F2bFe,F3bFe]

def F2Fe(x,Q2):
    return FFe(x,Q2)[0]

def F3Fe(x,Q2):
    return FFe(x,Q2)[1]

def F2bFe(x,Q2):
    return FFe(x,Q2)[2]

def F3bFe(x,Q2):
    return FFe(x,Q2)[3]

#Calculating the structure functions F2 and F3#

Q2a=np.linspace(1.69,136,1000) #limited kinematical space from EPPS21#
x=np.array([0.275,0.45]) #Insert desired value(s) of x here##
color=['crimson','steelblue'] 
for i in range (0,len(x)):
    T=F2Fe(x[i],Q2a)
    Tb=F2bFe(x[i],Q2a)

    S=F3Fe(x[i],Q2a)*x[i]
    Sb=F3bFe(x[i],Q2a)*x[i]
    plt.plot(Q2a,(T[0][:]+Tb[0][:])/2,color=color[i])
    plt.plot(Q2a,(S[0][:]+Sb[0][:])/2,color=color[i])