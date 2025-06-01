from parton import mkPDF
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import integrate as int


#Useful quantities-working in GeV#
pi=np.pi
Mp=9.3824*10**(-1) #proton mass#
Mtau=1.7769 #tau lepton mass#
Mmu=105.6584*10**(-3) #muon lepton mass#
MO=16*Mp #Oxygen mass#
MFe=56*Mp #Iron mass#
MW=80.3692 #W boson mass#
GF=1.1664*10**(-5) #Fermi's Constant#



##### PDFs--Structure Functions #####
pdfP = mkPDF('CT18ANLO',0) #making proton PDFs#
pdfO = mkPDF('EPPS21nlo_CT18Anlo_O16',0) #making Oxygen PDFs#
pdfFe = mkPDF('EPPS21nlo_CT18Anlo_Fe56',0) #making Iron PDFs#

#making proton structure functions for neutrinos#
def Fp(x,Q2):
    dp=pdfP.xfxQ2(1,x,Q2)
    dbp=pdfP.xfxQ2(-1,x,Q2)
    
    up=pdfP.xfxQ2(2,x,Q2)
    ubp=pdfP.xfxQ2(-2,x,Q2)
    
    sp=pdfP.xfxQ2(3,x,Q2)
    sbp=pdfP.xfxQ2(-3,x,Q2)
    
    cp=pdfP.xfxQ2(4,x,Q2)
    cbp=pdfP.xfxQ2(-4,x,Q2)
    
    bp=pdfP.xfxQ2(5,x,Q2)
    bbp=pdfP.xfxQ2(-5,x,Q2)
    
    F2p = 2*(dp + sp + ubp + bp + cbp)
    F2bp= 2*(up + cp + dbp + sbp + bbp)

    F3p = (2*(dp + sp - ubp + bp - cbp))/x
    F3bp= (2*(up + cp - dbp - sbp - bbp))/x 
    
    return [F2p,F3p,F2bp,F3bp]

def F2p(x,Q2):
    return Fp(x,Q2)[0]

def F3p(x,Q2):
    return Fp(x,Q2)[1]

def F2bp(x,Q2):
    return Fp(x,Q2)[2]

def F3bp(x,Q2):
    return Fp(x,Q2)[3]


#making Oxygen structure functions for neutrinos#
def FO(x,Q2):
    dO=pdfO.xfxQ2(1,x,Q2)
    dbO=pdfO.xfxQ2(-1,x,Q2)
    
    uO=pdfO.xfxQ2(2,x,Q2)
    ubO=pdfO.xfxQ2(-2,x,Q2)
    
    sO=pdfO.xfxQ2(3,x,Q2)
    sbO=pdfO.xfxQ2(-3,x,Q2)
    
    cO=pdfO.xfxQ2(4,x,Q2)
    cbO=pdfO.xfxQ2(-4,x,Q2)
    
    bO=pdfO.xfxQ2(5,x,Q2)
    bbO=pdfO.xfxQ2(-5,x,Q2)
    
    F2O = 2*(dO + sO + ubO + bO + cbO)
    F2bO= 2*(uO + cO + dbO + sbO + bbO)

    F3O = (2*(dO + sO - ubO + bO - cbO))/x
    F3bO= (2*(uO + cO - dbO - sbO - bbO))/x 
    
    return [F2O,F3O,F2bO,F3bO]

def F2O(x,Q2):
    return FO(x,Q2)[0]

def F3O(x,Q2):
    return FO(x,Q2)[1]

def F2bO(x,Q2):
    return FO(x,Q2)[2]

def F3bO(x,Q2):
    return FO(x,Q2)[3]



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



##### DIS Cross-Sections #####

def W(y,x,F2,F3,m,M,E):
    Q2=2*x*y*M*E
    
    W1 = ( x*y**2 + (m**2 * y)/(2*M*E) )*F2(x,Q2)/(2*x)
    W2 = ( 1 - y - (x*M*y)/(2*E) - (m**2)/(4*E**2) )*F2(x,Q2)
    W3 = ( x * y -(x * y**2)/2 - (m**2 * y)/(4*M*E) )*F3(x,Q2)
    W5 = ( -(m**2)/(2*M*E))*F2(x,Q2)/x
    
    return [W1,W2,W3,W5]

def W1(y,x,F2,F3,m,M,E):
    return W(y,x,F2,F3,m,M,E)[0]

def W2(y,x,F2,F3,m,M,E):
    return W(y,x,F2,F3,m,M,E)[1]

def W3(y,x,F2,F3,m,M,E):
    return W(y,x,F2,F3,m,M,E)[2]

def W5(y,x,F2,F3,m,M,E):
    return W(y,x,F2,F3,m,M,E)[3]

#neutrino diff cross-section#
def d2sdxdy(y,x,F2,F3,m,M,E):
    Q2=2*x*y*M*E
    d2sdxdy=((GF**2 * M * E)/pi) * (MW**4/((Q2 + MW**2)**2))*(W1(y,x,F2,F3,m,M,E) + W2(y,x,F2,F3,m,M,E) + W3(y,x,F2,F3,m,M,E) + W5(y,x,F2,F3,m,M,E))
    d2sdxdy=d2sdxdy*3.8939*10**(-28) #convert GeV-2 to cm2#
    return d2sdxdy

#anti-neutrino diff cross-section#
def d2bsdxdy(y,x,F2,F3,m,M,E):
    Q2=2*x*y*M*E
    d2bsdxdy=((GF**2 * M * E)/pi) * (MW**4/((Q2 + MW**2)**2))*(W1(y,x,F2,F3,m,M,E) + W2(y,x,F2,F3,m,M,E) - W3(y,x,F2,F3,m,M,E) + W5(y,x,F2,F3,m,M,E))
    d2bsdxdy=d2bsdxdy*3.8939*10**(-28) #convert GeV-2 to cm2#
    return d2bsdxdy


#neutrino cross-section#
def sE(F2,F3,m,M,E):
    Q2min=1.69
    xmin=1.69/(2*M*E)     ##These limits are fixed from the constraints of EPPS21##
    def ymin(x):
        ymin=1.69/(2*x*M*E)
        return ymin
    r=int.dblquad(d2sdxdy,xmin,1,ymin,1,args=(F2,F3,m,M,E))
    sE=r[0]
    return sE

#anti-neutrino cross-section#
def sbE(F2,F3,m,M,E):
    Q2min=1.69
    xmin=1.69/(2*M*E)
    def ymin(x):
        ymin=1.69/(2*x*M*E)
        return ymin
    r=int.dblquad(d2bsdxdy,xmin,1,ymin,1,args=(F2,F3,m,M,E))
    sbE=r[0]
    return sbE


##Considering H2O targets at relatively High Energy tau neutrinos##

E=np.logspace(4,np.log10(1.69*10**(7)/(2*Mp)),180) #The upper limit is fixed by the low bound x>=1e-7 of EPPS21#


apt=np.array([]) #proton-tau neutrino#
bpt=np.array([]) #proton-tau antineutrino#

aOt=np.array([]) #Oxygen-tau neutrino#
bOt=np.array([]) #Oxygen-tau antineutrino#


for i in range(0,len(E)):
    Spt=sE(F2p,F3p,Mtau,Mp,E[i])
    apt=np.append(apt,Spt)

    Sbpt=sbE(F2bp,F3bp,Mtau,Mp,E[i])
    bpt=np.append(bpt,Sbpt)

    SOt=sE(F2O,F3O,Mtau,Mp,E[i])
    aOt=np.append(aOt,SOt)

    SbOt=sbE(F2bO,F3bO,Mtau,Mp,E[i])
    bOt=np.append(bOt,SbOt)
    

#Building the neutrino-antineutrino average cross-section#

aH2Ot=(2*apt +16*aOt)/18
bH2Ot=(2*bpt +16*bOt)/18
                         #For tau neutrinos#
AH2Ot=(aH2Ot+bH2Ot)/2


### Plot the Results ###
plt.plot(E,AH2Ot,color='crimson',label=r'$\nu_{\tau}$')
plt.yscale("log")
plt.xscale("log")
plt.xlim(E[0],E[-1])
plt.minorticks_on()
plt.tick_params(direction='in',right='True',top='True',which='both')
plt.xlabel(r'$E$ [GeV]',fontsize=14)
plt.ylabel(r'$(\sigma^{\nu}_{H_2O} + \sigma^{\bar{\nu}}_{H_2O})/2 \, \, [ \text{cm}]$',fontsize=14)

##Considering Fe targets at relatively Low Energy muon neutrinos##
Eg=np.linspace(30,350,100)
nu=np.array([]) #neutrinos
bnu=np.array([]) #antineutrinos
for i in range(0,len(E)):
    nu=np.append(nu,sE(F2Fe,F3Fe,Mmu,Mp,Eg[i])*10**(38)/Eg[i])
    bnu=np.append(bnu,sbE(F2bFe,F3bFe,Mmu,Mp,Eg[i])*10**(38)/Eg[i]) #the crross-sections are normalised to the incoming neutrino Energy#

#Plot the results#
plt.plot(Eg,nu,color='crimson')
plt.plot(Eg,bnu,color='teal')

plt.minorticks_on()
plt.tick_params(direction='in',right='True',top='True',which='both')
plt.ylim(0.2,0.8)
plt.xlim(32,350)
plt.xlabel(r'$E$ [GeV]',fontsize=14)
plt.ylabel(r'$\sigma^{\nu}/E_{\nu} \, \, [10^{-38} \text{cm}^2/GeV]$',fontsize=14)
plt.annotate(r'$\nu$',xy=(80,0.60),color='crimson',fontsize=18)
plt.annotate(r'$\bar{\nu}$',xy=(80,0.23),color='teal',fontsize=18)