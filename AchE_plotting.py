import neuron
import matplotlib.pyplot as plt
from numpy import genfromtxt
import numpy as np
my_data = genfromtxt('Acetyl20Hz.csv', delimiter=',')
my_data=np.transpose(my_data)
print(my_data)
neuron.h.load_file('stdrun.hoc')
soma=neuron.h.Section(name='soma')
soma.insert('hh')
acetyl=neuron.h.concACh(soma(0.5))
ac_record=neuron.h.Vector()
ac_record.record(acetyl._ref_concentration)

ChIN_activity= np.concatenate([np.arange(0,20000,400),np.arange(20000,21000,50),np.arange(23000,30000,400),np.arange(30000,31000,50),np.arange(33000,40000,400),np.arange(40000,41000,50),np.arange(43000,50000,400),np.arange(50000,51000,50),np.arange(53000,60000,400),np.arange(60000,61000,50),np.arange(63000,70000,400),np.arange(70000,71000,50),np.arange(73000,80000,400)])
VecStim_ChIN=neuron.h.VecStim()
ChIN_Vector_Activity=neuron.h.Vector(ChIN_activity)
VecStim_ChIN.play(ChIN_Vector_Activity)


NETCON1=neuron.h.NetCon(VecStim_ChIN,acetyl)
NETCON1.delay=0
NETCON1.threshold=0
NETCON1.weight[0]=1

tSave=neuron.h.Vector()
tSave.record(neuron.h._ref_t)


neuron.h.tstop=50000
neuron.h.run()
import pdb
pdb.set_trace()
plt.plot(my_data[0][40:350]-1,my_data[1][40:350],'brown',label='M4 model')

plt.plot((np.array(tSave)[160000:1600000]-1150)*1e-3,(np.array(ac_record)[160000:1600000])*1e6,'orange',label="Model",linewidth=1.3)
plt.ylabel("Acetylcholine concentration (nM)")
plt.xlabel("Time (s)")

plt.legend()
plt.show()
plt.savefig('Output/Fig-5-M4/Mech-ACh.svg',format='svg')

