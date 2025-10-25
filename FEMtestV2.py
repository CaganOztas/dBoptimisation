import femm
import matplotlib.pyplot as plt
import time
import numpy as np
import scipy.ndimage, scipy.signal
femm.openfemm()
femm.newdocument(0)
femm.mi_probdef(0, 'millimeters', 'axi', 1.e-8, 0, 30)
length = 4
width = 4
innerradius = 1
Lsegments = 51
Bl = []

def generateCoil(width,length,innerRadius,wireAWG):
    femm.mi_drawrectangle(innerRadius,length/2,innerRadius+width,-length/2)
    femm.mi_addblocklabel(innerRadius+width/2 , 0)
    femm.mi_getmaterial(str(wireAWG)+" AWG")
    femm.mi_selectlabel(innerRadius+width/2 , 0)
    femm.mi_addcircprop("Coil", 1, 1)
    femm.mi_setblockprop(str(wireAWG)+" AWG", 0, 1, "Coil", 0, 0, 400)
    femm.mi_clearselected()
    
    return (innerRadius+width/2 , 0)

coilLabelPos = generateCoil(width,length,innerradius,22)
femm.mi_getmaterial("Air")
femm.mi_addblocklabel(0,0)
femm.mi_selectlabel(0,0)
femm.mi_setblockprop("Air")

femm.mi_makeABC()
femm.mi_saveas("TEST.fem")
femm.mi_analyze()
femm.mi_loadsolution()


def getCrossSectionB(z,segments):
    Bw = []
    dw = innerradius/(segments-1)
    integral = 0
    for i in range(segments):
        pos = innerradius - i * innerradius/(segments-1)
        B = femm.mo_getb(pos,z)[1]
        Bw.append(B)
        if(i > 0):
            segmentArea = (B+Bw[i-1])*dw/2 #
            integral += segmentArea
    return integral


def differentiate(values):
    delta = len
    rates = []
    for i in range(0,len(values)):
        if(i!=0 and i!=len(values)-1):
            rates.append( (values[i+1]-values[i-1])/(2*(length/(Lsegments-1))) )
        elif(i==0):
            rates.append((values[i+1]-values[i])/(length/(Lsegments-1)))
        else:
            rates.append((values[i]-values[i-1])/(length/(Lsegments-1)))
        print(length/(Lsegments-1))
    return rates
        
length += 10
positions = []
for i in range(Lsegments):
    pos = -length/2 + i * length/(Lsegments-1)
    positions.append(pos)
    B = getCrossSectionB(pos,15)
    Bl.append(B)
Bl_smoothed =  scipy.signal.savgol_filter(Bl, window_length=9, polyorder=3)

#gradient = differentiate(Bl_smoothed)
gradient = scipy.signal.savgol_filter(Bl, window_length=9, polyorder=3, deriv=1, delta=length/(Lsegments-1))
plt.scatter(positions, Bl, label = "B field", color='red', s=0.5)
plt.plot(positions, Bl_smoothed, label = "B field smoothed")
plt.plot(positions, gradient, label = "dB/dz")
plt.legend()
plt.plot(positions, [0 for i in range(len(positions))], 'k--')
plt.show()
#while True:
#    pass    