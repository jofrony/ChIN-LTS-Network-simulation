import elephant
import neuron
import neo
import quantities as pq
import matplotlib.pyplot as plt
from numpy import genfromtxt
import numpy as np

trace= np.loadtxt("Model_ChIN-Opt0-v9-voltage-cortex-10-thalamus-4-checker-0.03{'ChIN_TH': 140, 'LTS_Crtx': 66, 'ChIN_Crtx': 55}ChIN.txt")
traceLTS = np.loadtxt("Model_ChIN-Opt0-v9-voltage-cortex-10-thalamus-4-checker-0.03{'ChIN_TH': 140, 'LTS_Crtx': 66, 'ChIN_Crtx': 55}LTS.txt")

plt.figure(0)
plt.plot((trace[0][300000:540000]-7500)*1e-3,trace[1][300000:540000],label='Model v10',color='black')

plt.ylabel('Voltage (mV)')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig('SliceCHIN.svg',format='svg')
plt.show()
plt.clf()

plt.figure(1)
plt.plot((traceLTS[0][300000:540000]-7500)*1e-3,traceLTS[1][300000:540000],label='Model v10',color='black')

plt.ylabel('Voltage (mV)')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig('SliceLTS.svg',format='svg')
plt.show()

