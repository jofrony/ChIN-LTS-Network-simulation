import glob
import json
from Neuron_model_extended import NeuronModel
import sys
from neuron import h
import bluepyopt.ephys as ephys
import matplotlib.pyplot as plt
import numpy as np
import quantities as pq
import random
import time


with open("config-Figure6-cortical-synapse-ChIN-spiking.json",'r') as parameter_file:
        parameters_file = json.load(parameter_file)


start_time=time.time()

brainarea=parameters_file["Synapse"]["brainarea"]

synapsenumber="synapse_"+parameters_file["Synapse"]["synapse-parameter"][0]

with open("Synapses/cortical/cortical-synapse-model-parameter-ChIN.json",'r') as synapse_file:
    ChINsynapses =json.load(synapse_file)

synapse_ChIN=ChINsynapses[brainarea][synapsenumber]

synapsename=ChINsynapses[brainarea][synapsenumber]["type"]+"_"+ChINsynapses[brainarea][synapsenumber]["experiment"]

synapse_ChIN_conductance=parameters_file["Synapse"]["synapse-parameter"][1]


print("ChIN")

print(synapse_ChIN)

print(synapse_ChIN_conductance)


Cells = dict()

corticalsynapses = list()

configs= parameters_file["single-cell-models"] ##############

cellssynapseattached = ["Unknown[0]"]

position_cortical_input_ChIN=np.random.choice(range(397), 55, replace=False)

new_sim=ephys.simulators.NrnSimulator(cvode_active=False)

Cortical_activity=parameters_file["Cortical_stimulation"]

VecStim_Cortex=new_sim.neuron.h.VecStim()

Cortical_Vector_Activity=new_sim.neuron.h.Vector(Cortical_activity)

VecStim_Cortex.play(Cortical_Vector_Activity)

Gaussian = []

Synapse_TOTAL={"ChIN_TH":0,"LTS_Crtx":0,"ChIN_Crtx":0}

for cells in range(1):

    with open(configs[cells],'r') as config_file:
        config = json.load(config_file)

    morph = config["morphology"]
    print("morphology-DONE")

    param=config["parameters"]
    print("parameters-DONE")

    mech=config["mechanisms"]
    print("mechanisms done")

    info_cell = dict()

    ChIN_cell=NeuronModel(param_file=param,morph_file=morph,mech_file=mech)

    ChIN_cell.instantiate(sim=new_sim)

    info_cell.update({"soma_sec": ChIN_cell})
    
    info_cell.update({"soma_sim": new_sim})

    vSave = new_sim.neuron.h.Vector()

    for isec, sec in enumerate(ChIN_cell.icell.soma):

        for seg in sec:
            if "0.5" in str(seg):
                gauss=new_sim.neuron.h.InGauss(seg)
                gauss.delay=0
                gauss.dur=10e10
                gauss.mean=0
                gauss.stdev=0.05

                vSave.record(getattr(seg,'_ref_v'))
                info_cell.update({"soma_access":seg})
                spike_time = new_sim.neuron.h.Vector()
                recording_netcon = new_sim.neuron.h.NetCon(getattr(seg,'_ref_v'),None, sec = sec)
                recording_netcon.threshold = 0
                recording_netcon.record(spike_time)

    random_stream_Crtx=0

    for sec in new_sim.neuron.h.allsec():

        if cellssynapseattached[cells] in sec.name() and "soma" not in sec.name() and "axon" not in sec.name():

            for seg in sec:
                
                if random_stream_Crtx in position_cortical_input_ChIN:

                    stim_cortical=new_sim.neuron.h.tmGlut(seg)
                    stim_cortical.U=synapse_ChIN['U']
                    stim_cortical.tau=synapse_ChIN['tau']
                    stim_cortical.tauR=synapse_ChIN['tauR']
                    stim_cortical.tauR=synapse_ChIN['tauF']
                    stim_cortical.nmda_ratio=synapse_ChIN['nmda_ratio']

                    nc_cortical = new_sim.neuron.h.NetCon(VecStim_Cortex,stim_cortical)
                    nc_cortical.delay=1
                    nc_cortical.threshold=0
                    nc_cortical.weight[0]=synapse_ChIN_conductance

                    corticalsynapses.append([stim_cortical,nc_cortical])

                    Synapse_TOTAL['ChIN_Crtx']=Synapse_TOTAL['ChIN_Crtx']+1

                random_stream_Crtx=random_stream_Crtx+1

                
    info_cell.update({"spike_con": recording_netcon})
    info_cell.update({"spike_train": spike_time})

    info_cell.update({"soma_voltage":vSave})
    Cells.update({"ChIN_"+str(cells):info_cell})

print(Cells)

tSave = new_sim.neuron.h.Vector()
tSave.record(new_sim.neuron.h._ref_t)
print(Cortical_activity)

print(Synapse_TOTAL)
new_sim.neuron.h.tstop=parameters_file["simulation_time"]

new_sim.neuron.h.run()


k=0

name_cells={"ChIN_0":"ChIN"}

for cell,cell_info in Cells.items():
    print(cell)
    cell_v=cell_info["soma_voltage"]    
    plt.figure(k)
    #print(cell_info["vclamp"].i) 
    plt.plot(tSave,cell_v,label=name_cells[cell],c='black')
    np.savetxt("Output/Fig-6-Cortical/ChIN/Cortical-synapse-Opt0-v0s-spiking-"+brainarea+'-'+synapsenumber+'-'+synapsename+'-Number-of-synapses'+str(Synapse_TOTAL["ChIN_Crtx"])+'-'+name_cells[cell]+'.txt',[tSave,cell_v])
    plt.title("ChIN_LTS")
    plt.legend()
    plt.savefig("Output/Fig-6-Cortical/ChIN/Cortical-synapse-Opt0-v0s-spiking-"+brainarea+'-'+synapsenumber+'-'+synapsename+'-Number-of-synapses'+str(Synapse_TOTAL["ChIN_Crtx"])+'-'+name_cells[cell]+'.svg')
    plt.clf()
    k=k+1


print('wall time: {}s'.format(time.time() - start_time))





