import elephant
import neuron
import neo
import quantities as pq
import matplotlib.pyplot as plt
from numpy import genfromtxt
import numpy as np

Muscarinic_experiment_data = np.transpose(genfromtxt('../../Experimental-data/Muscarinic-effect-Melendez-Haidi.csv', delimiter=','))

v0_LTS= np.loadtxt('Model_LTS-M4-Model-LTS-v0.txt')
v1_LTS= np.loadtxt('Model_LTS-M4-Model-LTS-v1.txt')
v2_LTS= np.loadtxt('Model_LTS-M4-Model-LTS-v2.txt')
v3_LTS= np.loadtxt('Model_LTS-M4-Model-LTS-v3.txt')
v4_LTS= np.loadtxt('Model_LTS-M4-Model-LTS-v4.txt')
v5_LTS= np.loadtxt('Model_LTS-M4-Model-LTS-v5.txt')
v6_LTS= np.loadtxt('Model_LTS-M4-Model-LTS-v6.txt')
v7_LTS= np.loadtxt('Model_LTS-M4-Model-LTS-v7.txt')
v8_LTS= np.loadtxt('Model_LTS-M4-Model-LTS-v8.txt')
v9_LTS= np.loadtxt('Model_LTS-M4-Model-LTS-v9.txt')



plt.figure(0)
#plt.plot(v0_LTS[0],v0_LTS[1],label='Model v1',color='black')
#plt.plot(v1_LTS[0],v1_LTS[1],label='Model v2',color='black')
#plt.plot(v2_LTS[0],v2_LTS[1],label='Model v3',color='black')
#plt.plot(v3_LTS[0],v3_LTS[1],label='Model v4',color='black')
#plt.plot(v4_LTS[0],v4_LTS[1],label='Model v5',color='black')
#plt.plot(v5_LTS[0],v5_LTS[1],label='Model v6',color='black')
#plt.plot(v6_LTS[0],v6_LTS[1],label='Model v7',color='black')
#plt.plot(v7_LTS[0],v7_LTS[1],label='Model v8',color='black')
#plt.plot(v8_LTS[0],v8_LTS[1],label='Model v9',color='black')
plt.plot(v9_LTS[0][160000:320000]*1e-3,v9_LTS[1][160000:320000],label='Model v10',color='black')
plt.ylabel('Voltage (mV)')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig('Spiking-M4-LTS-v9.svg',format='svg')
plt.show()
plt.clf()

sts=list()

neo_voltage_v0=neo.AnalogSignal(np.array(v0_LTS[1]),sampling_period=0.025*pq.ms,units='mV')
spike_train_v0=elephant.spike_train_generation.peak_detection(neo_voltage_v0,threshold=-20 * pq.mV)
sts.append(neo.SpikeTrain(spike_train_v0,t_start=0 * pq.ms,t_stop=15000 * pq.ms))

neo_voltage_v1=neo.AnalogSignal(np.array(v1_LTS[1]),sampling_period=0.025*pq.ms,units='mV')
spike_train_v1=elephant.spike_train_generation.peak_detection(neo_voltage_v1,threshold=-20 * pq.mV)
sts.append(neo.SpikeTrain(spike_train_v1,t_start=0 * pq.ms,t_stop=15000 * pq.ms))

neo_voltage_v2=neo.AnalogSignal(np.array(v2_LTS[1]),sampling_period=0.025*pq.ms,units='mV')
spike_train_v2=elephant.spike_train_generation.peak_detection(neo_voltage_v2,threshold=-20 * pq.mV)
sts.append(neo.SpikeTrain(spike_train_v2,t_start=0 * pq.ms,t_stop=15000 * pq.ms))

neo_voltage_v3=neo.AnalogSignal(np.array(v3_LTS[1]),sampling_period=0.025*pq.ms,units='mV')
spike_train_v3=elephant.spike_train_generation.peak_detection(neo_voltage_v3,threshold=-20 * pq.mV)
sts.append(neo.SpikeTrain(spike_train_v3,t_start=0 * pq.ms,t_stop=15000 * pq.ms))

neo_voltage_v4=neo.AnalogSignal(np.array(v4_LTS[1]),sampling_period=0.025*pq.ms,units='mV')
spike_train_v4=elephant.spike_train_generation.peak_detection(neo_voltage_v4,threshold=-20 * pq.mV)
sts.append(neo.SpikeTrain(spike_train_v4,t_start=0 * pq.ms,t_stop=15000 * pq.ms))

neo_voltage_v5=neo.AnalogSignal(np.array(v5_LTS[1]),sampling_period=0.025*pq.ms,units='mV')
spike_train_v5=elephant.spike_train_generation.peak_detection(neo_voltage_v5,threshold=-20 * pq.mV)
sts.append(neo.SpikeTrain(spike_train_v5,t_start=0 * pq.ms,t_stop=15000 * pq.ms))

neo_voltage_v6=neo.AnalogSignal(np.array(v6_LTS[1]),sampling_period=0.025*pq.ms,units='mV')
spike_train_v6=elephant.spike_train_generation.peak_detection(neo_voltage_v6,threshold=-20 * pq.mV)
sts.append(neo.SpikeTrain(spike_train_v6,t_start=0 * pq.ms,t_stop=15000 * pq.ms))

neo_voltage_v7=neo.AnalogSignal(np.array(v7_LTS[1]),sampling_period=0.025*pq.ms,units='mV')
spike_train_v7=elephant.spike_train_generation.peak_detection(neo_voltage_v7,threshold=-20 * pq.mV)
sts.append(neo.SpikeTrain(spike_train_v7,t_start=0 * pq.ms,t_stop=15000 * pq.ms))

neo_voltage_v8=neo.AnalogSignal(np.array(v8_LTS[1]),sampling_period=0.025*pq.ms,units='mV')
spike_train_v8=elephant.spike_train_generation.peak_detection(neo_voltage_v8,threshold=-20 * pq.mV)
sts.append(neo.SpikeTrain(spike_train_v8,t_start=0 * pq.ms,t_stop=15000 * pq.ms))

neo_voltage_v9=neo.AnalogSignal(np.array(v9_LTS[1]),sampling_period=0.025*pq.ms,units='mV')
spike_train_v9=elephant.spike_train_generation.peak_detection(neo_voltage_v9,threshold=-20 * pq.mV)
sts.append(neo.SpikeTrain(spike_train_v9,t_start=0 * pq.ms,t_stop=15000 * pq.ms))

binsize = 250 * pq.ms
pop_count = elephant.statistics.time_histogram(sts, binsize)

plt.figure(1)

mean=np.average(np.transpose(pop_count)[0][0:60])

plt.figure(0)
plt.bar(Muscarinic_experiment_data[0],Muscarinic_experiment_data[1],width=0.25,label='Experimental data',fill=False,edgecolor='red')
plt.xlabel('Time (s)')
plt.ylabel('Normalized spike rate')
plt.legend()
plt.savefig('Spiking-Normalized-M4-LTS-experiment.svg',format='svg')

plt.bar(np.array(Muscarinic_experiment_data[0]),np.transpose(pop_count)[0]*(1/mean),label='Model',fill=False,edgecolor='black',width=0.25)
plt.xlabel('Time (s)')
plt.ylabel('Normalized spike rate')
plt.legend()
plt.savefig('Spiking-Normalized-M4-LTS.svg',format='svg')
plt.show()
