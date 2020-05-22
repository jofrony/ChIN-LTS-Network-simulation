import json
from Neuron_model_extended import NeuronModel
import sys
from neuron import h
import bluepyopt.ephys as ephys
import matplotlib.pyplot as plt
import numpy as np
import random
import time


with open("config-Figure6-thalamic-synapse-ChIN.json",'r') as parameter_file:
        parameters_file = json.load(parameter_file)

start_time=time.time()


brainarea=parameters_file["Synapse"]["brainarea"]

synapsenumber="synapse_"+parameters_file["Synapse"]["synapse-parameter"][0]

with open("Synapses/thalamic/thalamic-synapse-model-parameter-ChIN.json",'r') as synapse_file:
    ChINsynapses =json.load(synapse_file)


synapse_ChIN=ChINsynapses[brainarea][synapsenumber]

synapsename=ChINsynapses[brainarea][synapsenumber]["type"]+"_"+ChINsynapses[brainarea][synapsenumber]["experiment"]

synapse_ChIN_conductance=parameters_file["Synapse"]["synapse-parameter"][1]

position_for_thalamic_input=np.random.choice(range(395), 140, replace=False)

holding=parameters_file["Holding_voltage"]


print("ChIN")

print(synapse_ChIN)

print(synapse_ChIN_conductance)

Cells = dict()

thalamicsynapses = list()

configs= parameters_file["single-cell-models"] ##############

cellssynapseattached = ["Unknown[0]"]

new_sim=ephys.simulators.NrnSimulator(cvode_active=False)

Thalamic_activity= parameters_file["Thalamic_stimulation"]

VecStim_TH=new_sim.neuron.h.VecStim()

Thalamic_Vector_Activity=new_sim.neuron.h.Vector(Thalamic_activity)

VecStim_TH.play(Thalamic_Vector_Activity)

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
    print("print mechanisms done")

    info_cell = dict()

    ChIN_cell=NeuronModel(param_file=param,morph_file=morph,mech_file=mech)
    new_sim=ephys.simulators.NrnSimulator(cvode_active=False)
    ChIN_cell.instantiate(sim=new_sim)



    info_cell.update({"soma_sec": ChIN_cell})
    info_cell.update({"soma_sim": new_sim})

    vSave = new_sim.neuron.h.Vector()

    for isec, sec in enumerate(ChIN_cell.icell.soma):

        for seg in sec:
            if "0.5" in str(seg):

                stim_v=new_sim.neuron.h.IClamp(seg)
                stim_v.amp=holding
                stim_v.dur = 1e9
                stim_v.delay=0
                vSave.record(getattr(seg,'_ref_v'))
                info_cell.update({"soma_access":seg})
                info_cell.update({"vclamp":stim_v})
                spike_time = new_sim.neuron.h.Vector()
                recording_netcon = new_sim.neuron.h.NetCon(getattr(seg,'_ref_v'),None, sec = sec)
                recording_netcon.threshold = 0
                recording_netcon.record(spike_time)


    random_stream_TH=0

    for sec in new_sim.neuron.h.allsec():

        if cellssynapseattached[cells] in sec.name() and "axon" not in sec.name():

            for seg in sec:

                if random_stream_TH in position_for_thalamic_input:
                    
                    stim_thalamic=new_sim.neuron.h.tmGlut(seg)
                    stim_thalamic.U=synapse_ChIN["U"]
                    stim_thalamic.tau=synapse_ChIN['tau']
                    stim_thalamic.tauR=synapse_ChIN['tauR']
                    stim_thalamic.tauR=synapse_ChIN['tauF']
                    stim_thalamic.nmda_ratio=synapse_ChIN['nmda_ratio']

                    nc_thalamic = new_sim.neuron.h.NetCon(VecStim_TH,stim_thalamic)
                    nc_thalamic.delay=1
                    nc_thalamic.threshold=0
                    nc_thalamic.weight[0]=synapse_ChIN_conductance

                    thalamicsynapses.append([stim_thalamic,nc_thalamic])

                    Synapse_TOTAL['ChIN_TH']=Synapse_TOTAL['ChIN_TH']+1

                random_stream_TH=random_stream_TH+1

    info_cell.update({"spike_con": recording_netcon})
    
    info_cell.update({"spike_train": spike_time})

    info_cell.update({"soma_voltage":vSave})
    
    Cells.update({"ChIN_"+str(cells):info_cell})


tSave = new_sim.neuron.h.Vector()

tSave.record(new_sim.neuron.h._ref_t)

print(Thalamic_activity)

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
    np.savetxt("Output/Fig-6-Thalamic/Thalamic-synapse-"+parameters_file["single-cell-models"][0].split("/")[-1].split(".")[0]+brainarea+'-'+synapsenumber+'-'+synapsename+'-Number-of-synapses'+str(Synapse_TOTAL["ChIN_TH"])+'-'+name_cells[cell]+'.txt',[tSave,cell_v])
    plt.title("ChIN_LTS")
    plt.legend()
    plt.savefig("Output/Fig-6-Thalamic/Thalamic-synapse-"+parameters_file["single-cell-models"][0].split("/")[-1].split(".")[0]+brainarea+'-'+synapsenumber+'-'+synapsename+'-Number-of-synapses'+str(Synapse_TOTAL["ChIN_TH"])+'-'+name_cells[cell]+'.svg')
    plt.clf()
    
    k=k+1

print('wall time: {}s'.format(time.time() - start_time))

        


