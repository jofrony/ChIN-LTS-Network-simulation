import glob
import json
from Neuron_model_extended import NeuronModel
import sys
import numpy as np
from neuron import h
import bluepyopt.ephys as ephys 
import matplotlib.pyplot as plt
import random
import time

with open("config-Figure6-cortical-synapse-LTS-spiking.json",'r') as parameter_file:
        parameters_file = json.load(parameter_file)


start_time=time.time()

brainarea=parameters_file["Synapse"]["brainarea"]

synapsenumber="synapse_"+parameters_file["Synapse"]["synapse-parameter"][0]

with open("Synapses/cortical/cortical-synapse-model-parameter-LTS.json",'r') as synapse_file:
    LTSsynapses =json.load(synapse_file)

synapse_LTS=LTSsynapses[brainarea][synapsenumber]

synapsename=LTSsynapses[brainarea][synapsenumber]["type"]+"_"+LTSsynapses[brainarea][synapsenumber]["experiment"]

synapse_LTS_conductance=parameters_file["Synapse"]["synapse-parameter"][1]

position_for_cortical_input_LTS=np.arange(66)

print("LTS")

print(synapse_LTS)

print(synapse_LTS_conductance)

Cells = dict()

corticalsynapses = list()

configs= parameters_file["single-cell-models"] ######################

cellssynapseattached = ["Unknown[0]"]

new_sim=ephys.simulators.NrnSimulator(cvode_active=False)

Cortical_activity = parameters_file["Cortical_stimulation"]

VecStim_Cortex=new_sim.neuron.h.VecStim()

Cortical_Vector_Activity=new_sim.neuron.h.Vector(Cortical_activity)

VecStim_Cortex.play(Cortical_Vector_Activity)

Gaussian = []

Synapse_TOTAL={"LTS_Crtx":0}

for cells in range(1):

    with open(configs[cells],'r') as config_file:
        config = json.load(config_file)

    morph = config["morphology"]
    print("morphology-DONE")

    param=config["parameters"]
    print("parameters-DONE")

    mech=config["mechanisms"]
    print("print mechanisms done")

    info_cell = dict()

    LTS_cell=NeuronModel(param_file=param,morph_file=morph,mech_file=mech)

    LTS_cell.instantiate(sim=new_sim)

    info_cell.update({"soma_sec": LTS_cell})
    info_cell.update({"soma_sim": new_sim})

    vSave = new_sim.neuron.h.Vector()

    for isec, sec in enumerate(LTS_cell.icell.soma):

        for seg in sec:
            if "0.5" in str(seg):
                gauss=new_sim.neuron.h.InGauss(seg)
                gauss.delay=0
                gauss.dur=10e10
                gauss.mean=0
                gauss.stdev=0.05
                Gaussian.append(gauss)
          
                vSave.record(getattr(seg,'_ref_v'))
                info_cell.update({"soma_access":seg})
                
                spike_time = new_sim.neuron.h.Vector()
                recording_netcon = new_sim.neuron.h.NetCon(getattr(seg,'_ref_v'),None, sec = sec)
                recording_netcon.threshold = 0
                recording_netcon.record(spike_time)

    random_stream_Crtx_LTS=0
    
    for sec in new_sim.neuron.h.allsec():

        if cellssynapseattached[cells] in sec.name() and "axon" not in sec.name():

            for seg in sec:

                for k in range(4):

                        
                        if random_stream_Crtx_LTS in position_for_cortical_input_LTS:
                
                            #print('attaching LTS Cortical')
                            #print(cellssynapseattached[cells])
                            stim_cortical=new_sim.neuron.h.tmGlut(0.5,sec=sec)
                            stim_cortical.U=synapse_LTS['U']
                            stim_cortical.tau=synapse_LTS['tau']
                            stim_cortical.tauR=synapse_LTS['tauR']
                            stim_cortical.tauR=synapse_LTS['tauF']
                            stim_cortical.nmda_ratio=synapse_LTS['nmda_ratio']

                            nc_cortical = new_sim.neuron.h.NetCon(VecStim_Cortex,stim_cortical)
                            nc_cortical.delay=1
                            nc_cortical.threshold=0
                            nc_cortical.weight[0]=synapse_LTS_conductance

                            corticalsynapses.append([stim_cortical,nc_cortical])

                            Synapse_TOTAL['LTS_Crtx']=Synapse_TOTAL['LTS_Crtx']+1

                        random_stream_Crtx_LTS=random_stream_Crtx_LTS+1


                   

    info_cell.update({"spike_con": recording_netcon})
    
    info_cell.update({"spike_train": spike_time})

    info_cell.update({"soma_voltage":vSave})
    
    Cells.update({"LTS_"+str(cells):info_cell})
    
print(Cells)

allsynapses=list()

Muscarinic=True

Muscarinic_true = list()

for sec in new_sim.neuron.h.allsec():
    
    for segi in sec:
        
        if "axon" not in str(segi) and Muscarinic==True:
            
            print('Muscarinic intracellular')
            #print(segi)
            ach_conce = new_sim.neuron.h.concACh(segi)
            pointer_concentration=ach_conce._ref_concentration
            Minput = new_sim.neuron.h.M4(segi)
            new_sim.neuron.h.setpointer(pointer_concentration,'conc_ACH',Minput)
            pointer_modulation=Minput._ref_Ach_M4R
            new_sim.neuron.h.setpointer(pointer_modulation,'muscarinic_modulation', segi.kir23_lts)
            Muscarinic_true.append([ach_conce,Minput,pointer_modulation,pointer_concentration])
            muscarinic_recording= new_sim.neuron.h.Vector()
            muscarinic_recording.record(Minput._ref_Ach_M4R)
            acetyl_recording= new_sim.neuron.h.Vector()
            acetyl_recording.record(ach_conce._ref_concentration)


tSave = new_sim.neuron.h.Vector()

tSave.record(new_sim.neuron.h._ref_t)

print(Cortical_activity)

print(Synapse_TOTAL)

new_sim.neuron.h.tstop=parameters_file["simulation_time"]

new_sim.neuron.h.run()

k=0

name_cells={"LTS_0":"LTS"}

for cell,cell_info in Cells.items():
    print(cell)
    cell_v=cell_info["soma_voltage"]    
    plt.figure(k)
    #print(cell_info["vclamp"].i) 
    plt.plot(tSave,cell_v,label=name_cells[cell],c='black')
    np.savetxt("Output/Fig-6-Cortical/LTS/Cortical-synapse-spiking-"+brainarea+'-'+synapsenumber+'-'+synapsename+'-Number-of-synapses'+str(Synapse_TOTAL["LTS_Crtx"])+'-'+name_cells[cell]+'.txt',[tSave,cell_v])
    plt.title("ChIN_LTS")
    plt.legend()
    plt.savefig("Output/Fig-6-Cortical/LTS/Cortical-synapse-spiking-"+brainarea+'-'+synapsenumber+'-'+synapsename+'-Number-of-synapses'+str(Synapse_TOTAL["LTS_Crtx"])+'-'+name_cells[cell]+'.svg')
    plt.clf()
    k=k+1


print('wall time: {}s'.format(time.time() - start_time))





