import elephant
import neuron
import neo
import quantities as pq
import matplotlib.pyplot as plt
from numpy import genfromtxt
import numpy as np

traceLTS= np.loadtxt("Model_ChIN-Opt0-v9-voltage-cortex-20-thalamus-4-checker-0.05{'ChIN_TH': 140, 'LTS_Crtx': 66, 'ChIN_Crtx': 55}ChIN.txt")
traceChIN= np.loadtxt("Model_ChIN-Opt0-v9-voltage-cortex-20-thalamus-4-checker-0.05{'ChIN_TH': 140, 'LTS_Crtx': 66, 'ChIN_Crtx': 55}LTS.txt")

plt.plot((traceLTS[0][240000:1200000]-6000)*1e-3,traceLTS[1][240000:1200000],label='S5',color='black')
plt.ylabel('Voltage (mV)')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig('SliceLTS.svg',format='svg')
plt.show()
plt.clf()

plt.plot((traceChIN[0][240000:1200000]-6000)*1e-3,traceChIN[1][240000:1200000],label='S5',color='black')
plt.ylabel('Voltage (mV)')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig('SliceChIN.svg',format='svg')
plt.show()
