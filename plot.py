# Imports
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import *
import scipy as spi


# Reads in Data as 2x2 array
data = genfromtxt('Orion2ChemOut.txt',unpack=True)

#doesn't include last end point in array. If you want the first five
#items, do 0:6.
time = data[0,:]
HP = data[1,:]
H = data[2,:]
HM = data[3,:]
H2P = data[4,:]
H2 = data[5,:]
DP = data[6,:]
D = data[7,:]
DM = data[8,:]
HDP = data[9,:]
HD = data[10,:]
D2P = data[11,:]
D2 = data[12,:]
HEP = data[13,:]
HE = data[14,:]
HEPP = data[15,:]
ELEC = data[16,:]

Dens = data[17,:]
Temp = data[18,:]
EngD = data[19,:]
Pres = data[20,:]
mfrac = data[21,:]
gamma = data[22,:]

Pres = Pres/Pres[0]

time = log10(time)
#x = array([2,8,32,64,96,128,160,192])
#y = array([0.583,0.523,1.071,1.364,1.002,0.936,0.91,0.897])



#plt.setp(ax,xticklabels=[]) #turns off labels

# ---===---===---===---===---===---===---===---===---===
plt.subplot(1,1,1)
plt.subplots_adjust(left=0.18,bottom=0.15,wspace=0.0001,hspace=0.0001)

#plt.axes([0.16,0.13,0.79,0.82]) #left,bot,wid,height
plt.axis([1.0, 15, 1e-10, 1e6])
plt.axis([0.10, 13, 1e-20, 1e5])

# Axis labels
plt.xlabel(r'Log$_{10}$(Time)     (sec)') #TeX requires the r
plt.ylabel(r'Mass Fraction')

#Moves the x-axis down a little
plt.gca().xaxis.labelpad=9 #5 is default

#Makes ticks thicker/longer
ax = plt.gca()
ax.tick_params(which='both',width=2)
ax.tick_params(which='major',size=12)
ax.tick_params(which='minor',size=7)
	 
#moves numbers out a little
ax.tick_params(axis='y',pad=9)
ax.tick_params(axis='x',pad=5)


# Adds dashed lines
#xline = [-5,5]
#yline = [2,2]
#plt.plot(xline,yline,'k:',linewidth=2.0)

dred = [0.6,0,0]
dblue = [0,0,0.6]
dgreen = [0,0.6,0]
dpurp = [150.0/255.0,50.0/255.0,255.0/255.0]
dyel = [204.0/255.0,204.0/255.0,0/255.0]
# Plots data

plt.semilogy(time,HP,'--',color='black',label=r"H$^+$")
plt.semilogy(time,H,color='black',label=r"H")
plt.semilogy(time,HM,'-.',color='black',label=r"H$^-$")
plt.semilogy(time,H2P,'--',color=dyel,label=r"H$_2$$^+$")
plt.semilogy(time,H2,color=dyel,label=r"H$_2$")
plt.semilogy(time,DP,'--',color=dred,label=r"D$^+$")
plt.semilogy(time,D,color=dred,label=r"D")
plt.semilogy(time,DM,'-.',color=dred,label=r"D$^-$")
plt.semilogy(time,HDP,'--',color=dblue,label=r"HD$^+$")
plt.semilogy(time,HD,color=dblue,label=r"HD")
plt.semilogy(time,D2P,'--',color=dgreen,label=r"D$_2$$^+$")
plt.semilogy(time,D2,color=dgreen,label=r"D$_2$")
plt.semilogy(time,HEP,'--',color=dpurp,label=r"HE$^+$")
plt.semilogy(time,HE,color=dpurp,label=r"HE")
plt.semilogy(time,HEPP,'-.',color=dpurp,label=r"HE$^{++}$")
plt.semilogy(time,ELEC,':',color='black',label=r"e")

plt.semilogy(time,Temp,':',color=dred,label=r"Temp")
plt.semilogy(time,Pres,':',color=dpurp,label=r"Pres (n)")
plt.semilogy(time,gamma,':',color=dyel,label=r"Gamma")

#both with 1 particle / cm^3
analyticH2 = 1/(1+4*1e-13*10**(time)) #ionized to neutral H
analyticH2 = exp(-10**(time)*1300000000*1e-21) #photodissociate with J21=1


#plt.semilogy(time,analyticH2,':',color=dpurp,label=r"Analytic")


#plt.semilogx(array([164,164]),array([0.2,0.6]),'--',color='black')



#plt.semilogy(LineX,BH447*LineY,'--',color=dred)


#Text
#plt.text(1.4,1.03,r'$\beta = 1.0$ , ${\cal M}=1.41$',fontsize=24.0)

plt.legend(bbox_to_anchor=(0., 1.02, 1., .12),loc=1,numpoints=1,ncol=7,labelspacing=0.1,columnspacing=0.1,fontsize=10)

# Creates Plot
#plt.show()
plt.savefig('ChemOutput.pdf')
#plt.savefig('figure16.eps')

#Restores defaults
mpl.rcdefaults() 