from Neuron_model_extended import NeuronModel
import sys
from neuron import h
import bluepyopt.ephys as ephys
import json     
import matplotlib.pyplot as plt
import numpy as np
import elephant
import neo
import quantities as pq
import random
import time


with open("config-S5.json",'r') as parameter_file:
        parameters_file = json.load(parameter_file)


np.random.seed(2020)

# 2000 ms is 120 s

simulation_time= parameters_file["simulation_time"]


new_sim=ephys.simulators.NrnSimulator(dt=0.025,cvode_active=False)

cortex_frequency=parameters_file["Cortex_frequency"]
thalamus_frequency=parameters_file["Thalamus_frequency"]
checkerweight=parameters_file["checkerweight"]

Checkerweight=[1.1,checkerweight] #0.05 for LTS there is no response

Background_activity_Th = dict()

for i in range(752):
    np.random.seed(i)
    Th_B = elephant.spike_train_generation.homogeneous_poisson_process(4* pq.Hz,t_stop = simulation_time * pq.ms, t_start=0 * pq.ms)
    Vec_Th_B=new_sim.neuron.h.VecStim()
    Th_A_B=new_sim.neuron.h.Vector(Th_B)
    Vec_Th_B.play(Th_A_B)
    Background_activity_Th.update({"input_"+str(i): Vec_Th_B})
    
position_for_thalamic_input=np.random.choice(range(752), 140, replace=False)

Burst_activity_Th = dict()

for position in position_for_thalamic_input:
    Thalamic_activity = []
    VecStim_TH=new_sim.neuron.h.VecStim()
    Thalamic_Vector_Activity=new_sim.neuron.h.Vector(Thalamic_activity)
    VecStim_TH.play(Thalamic_Vector_Activity)
    Burst_activity_Th.update({"input_"+str(position): VecStim_TH})



'''

Cortical input 

'''

Background_activity_Crtx_ChIN = dict()

Background_activity_Crtx_LTS = dict()

position_background_for_cortical_input_LTS=np.arange(250)

position_background_for_cortical_input_ChIN=np.random.choice(range(752),294, replace=False)
                                                 
print(position_background_for_cortical_input_ChIN)
                                                 
for i in range(752):
    np.random.seed(i)
    Crtx_B = elephant.spike_train_generation.homogeneous_poisson_process(2* pq.Hz,t_stop = simulation_time * pq.ms, t_start=0 * pq.ms)
    Vec_Crtx_B=new_sim.neuron.h.VecStim()
    Crtx_A_B=new_sim.neuron.h.Vector(Crtx_B)
    Vec_Crtx_B.play(Crtx_A_B)

    if i in position_background_for_cortical_input_LTS:
        Background_activity_Crtx_LTS.update({"input_"+str(i): Vec_Crtx_B})
        
    if i in position_background_for_cortical_input_ChIN:
        Background_activity_Crtx_ChIN.update({"input_"+str(i): Vec_Crtx_B})
        

print(Background_activity_Crtx_ChIN)
print(Background_activity_Crtx_LTS)

position_cortical_input_ChIN=np.random.choice(range(395), 55, replace=False)
#position_cortical_input_ChIN=np.random.choice(range(397), 100, replace=False)

position_cortical_input_LTS=np.arange(66)

activity_Crtx_ChIN=dict()

activity_Crtx_LTS=dict()

Crtx_Stim = elephant.spike_train_generation.homogeneous_poisson_process(cortex_frequency* pq.Hz,t_stop = 6800 * pq.ms, t_start=6500 * pq.ms)

Vec_Crtx=new_sim.neuron.h.VecStim()
Crtx_A=new_sim.neuron.h.Vector(Crtx_Stim)
Vec_Crtx.play(Crtx_A)

for i in range(752):
    

    if i in position_cortical_input_ChIN:
        activity_Crtx_ChIN.update({"input_"+str(i): Vec_Crtx})

    if i in position_cortical_input_LTS:
        activity_Crtx_LTS.update({"input_"+str(i): Vec_Crtx})
        
    
start_time=time.time()

#From Doig et al 2014 8450, std=694 total amount of synapses
#752 from Thalamus
#294 from Cortex
#Inhibitory 5166

#ChIN index 0
#LTS index 1


Cells = dict()

NOsynapses = list()

thalamicsynapses = list()

corticalsynapses = list()

muscsynapses=list()


configs= parameters_file["single-cell-models"]

cellssynapseattached = ["Unknown[0]", "Unknown[1]"]

synapse_conductance_ChIN_cortical= parameters_file["Cortex_synapse_ChIN_conductance"]

synapse_conductance_ChIN_thalamic=parameters_file["Thalamus_synapse_ChIN_conductance"]

synapse_conductance_LTS_cortical = parameters_file["Cortex_synapse_LTS_conductance"]

synapse_config_cortical_MR = parameters_file["Cortex-synapse"]

synapse_config_thalamic = parameters_file["Thalamus-synapse"]

Gaussian = []

Synapse_TOTAL={"ChIN_TH":0,"LTS_Crtx":0,"ChIN_Crtx":0}

name= ["ChIN","LTS"]
print()


for cells in range(2):

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

    random_stream_TH=0
    random_stream_Crtx_ChIN=0
    random_stream_Crtx_LTS=0
    
    for sec in new_sim.neuron.h.allsec():
        

        if cellssynapseattached[cells] in sec.name() and "Unknown[0]" in sec.name() and "soma" not in sec.name() and "axon" not in sec.name():

            for seg in sec:

                for k in range(2):

                    if random_stream_TH in range(752):
                
                        print('attaching ChIN thalamic')
                        stim_thalamic_B=new_sim.neuron.h.tmGlut(seg)
                        synapse=synapse_config_thalamic["Unknown[0]"]
                        stim_thalamic_B.U=synapse["U"]
                        stim_thalamic_B.tau=synapse['tau']
                        stim_thalamic_B.tauR=synapse['tauR']
                        stim_thalamic_B.tauR=synapse['tauF']
                        stim_thalamic_B.nmda_ratio=synapse['nmda_ratio']

                        stimulation=Background_activity_Th["input_"+str(random_stream_TH)]
                        nc_thalamic_B = new_sim.neuron.h.NetCon(stimulation,stim_thalamic_B)
                        nc_thalamic_B.delay=1
                        nc_thalamic_B.threshold=0
                        nc_thalamic_B.weight[0]=synapse_conductance_ChIN_thalamic


                        thalamicsynapses.append([stim_thalamic_B,nc_thalamic_B])

                    if random_stream_TH in position_for_thalamic_input:


                        #print('attaching ChIN thalamic')
                        #print(cellssynapseattached[cells])
                        stim_thalamic=new_sim.neuron.h.tmGlut(seg)
                        synapse=synapse_config_thalamic["Unknown[0]"]
                        stim_thalamic.U=synapse["U"]
                        stim_thalamic.tau=synapse['tau']
                        stim_thalamic.tauR=synapse['tauR']
                        stim_thalamic.tauR=synapse['tauF']
                        stim_thalamic.nmda_ratio=synapse['nmda_ratio']

                        stimulation_burst=Burst_activity_Th["input_"+str(random_stream_TH)]

                        nc_thalamic = new_sim.neuron.h.NetCon(stimulation_burst,stim_thalamic)
                        nc_thalamic.delay=1
                        nc_thalamic.threshold=0
                        nc_thalamic.weight[0]=synapse_conductance_ChIN_thalamic

                        Synapse_TOTAL['ChIN_TH']=Synapse_TOTAL['ChIN_TH']+1

                        thalamicsynapses.append([stim_thalamic,nc_thalamic])

                    if random_stream_Crtx_ChIN in position_background_for_cortical_input_ChIN:


                        print('attaching ChIN Cortical')
                        #print(cellssynapseattached[cells])

                        stim_cortical_B=new_sim.neuron.h.tmGlut(seg)
                        synapse=synapse_config_cortical_MR["Unknown[0]"]
                        stim_cortical_B.U=synapse['U']
                        stim_cortical_B.tau=synapse['tau']
                        stim_cortical_B.tauR=synapse['tauR']
                        stim_cortical_B.tauR=synapse['tauF']
                        stim_cortical_B.nmda_ratio=synapse['nmda_ratio']


                        stimulation_Crtx_ChIN_B=Background_activity_Crtx_ChIN["input_"+str(random_stream_Crtx_ChIN)]

                        nc_cortical_B = new_sim.neuron.h.NetCon(stimulation_Crtx_ChIN_B,stim_cortical_B)
                        nc_cortical_B.delay=1
                        nc_cortical_B.threshold=0
                        nc_cortical_B.weight[0]=synapse_conductance_ChIN_cortical

                        corticalsynapses.append([stim_cortical_B,nc_cortical_B])

                    if random_stream_Crtx_ChIN in position_cortical_input_ChIN:



                        stim_cortical=new_sim.neuron.h.tmGlut(seg)
                        synapse=synapse_config_cortical_MR["Unknown[0]"]
                        stim_cortical.U=synapse['U']
                        stim_cortical.tau=synapse['tau']
                        stim_cortical.tauR=synapse['tauR']
                        stim_cortical.tauR=synapse['tauF']
                        stim_cortical.nmda_ratio=synapse['nmda_ratio']


                        stimulation_Crtx_ChIN=activity_Crtx_ChIN["input_"+str(random_stream_Crtx_ChIN)]

                        nc_cortical = new_sim.neuron.h.NetCon(stimulation_Crtx_ChIN,stim_cortical)
                        nc_cortical.delay=1
                        nc_cortical.threshold=0
                        nc_cortical.weight[0]=synapse_conductance_ChIN_cortical

                        corticalsynapses.append([stim_cortical,nc_cortical])

                        Synapse_TOTAL['ChIN_Crtx']=Synapse_TOTAL['ChIN_Crtx']+1
                        
                    random_stream_TH=random_stream_TH+1
                    random_stream_Crtx_ChIN=random_stream_Crtx_ChIN+1
              

        elif cellssynapseattached[cells] in sec.name() and "Unknown[1]" in sec.name() and "axon" not in sec.name():

            for seg in sec:

                for k in range(4):
                    
                    if random_stream_Crtx_LTS in position_background_for_cortical_input_LTS:
                        
                        stim_cortical_B=new_sim.neuron.h.tmGlut(seg)
                        synapse=synapse_config_cortical_MR["Unknown[1]"]
                        stim_cortical_B.U=synapse['U']
                        stim_cortical_B.tau=synapse['tau']
                        stim_cortical_B.tauR=synapse['tauR']
                        stim_cortical_B.tauR=synapse['tauF']
                        stim_cortical_B.nmda_ratio=synapse['nmda_ratio']

                        stimulation_B=Background_activity_Crtx_LTS["input_"+str(random_stream_Crtx_LTS)]
                
                        nc_cortical_B = new_sim.neuron.h.NetCon(stimulation_B,stim_cortical_B)
                        nc_cortical_B.delay=1
                        nc_cortical_B.threshold=0
                        nc_cortical_B.weight[0]=synapse_conductance_LTS_cortical

                        
                        corticalsynapses.append([stim_cortical_B,nc_cortical_B])

                    if random_stream_Crtx_LTS in position_cortical_input_LTS:


                        stim_cortical=new_sim.neuron.h.tmGlut(seg)
                        synapse=synapse_config_cortical_MR["Unknown[1]"]
                        stim_cortical.U=synapse['U']
                        stim_cortical.tau=synapse['tau']
                        stim_cortical.tauR=synapse['tauR']
                        stim_cortical.tauR=synapse['tauF']
                        stim_cortical.nmda_ratio=synapse['nmda_ratio']

                        stimulation=activity_Crtx_LTS["input_"+str(random_stream_Crtx_LTS)]

                        nc_cortical = new_sim.neuron.h.NetCon(stimulation,stim_cortical)
                        nc_cortical.delay=1
                        nc_cortical.threshold=0
                        nc_cortical.weight[0]=synapse_conductance_LTS_cortical


                        corticalsynapses.append([stim_cortical,nc_cortical])
                        Synapse_TOTAL['LTS_Crtx']=Synapse_TOTAL['LTS_Crtx']+1





                    random_stream_Crtx_LTS=random_stream_Crtx_LTS+1
                    


    info_cell.update({"spike_con": recording_netcon})
    info_cell.update({"spike_train": spike_time})

    info_cell.update({"soma_voltage":vSave})
    Cells.update({"ChIN_"+str(cells):info_cell})

'''
#######################################


Implementing the network



########################################
'''

allsynapses=list()


"CONSTRUCTING THE NETWORK"
Muscarinic=True
NO=True

for sec in new_sim.neuron.h.allsec():
    
    for seg in sec:

        if str(Cells['ChIN_0']['soma_access']) in str(seg):

            #print(str(Cells['ChIN_0']['soma_access']))
            print('ATTACHING CHIN SPIKING TO THE LTS')

            ChINprojectingsynapses = list()
            ChINnicotinicsynapses = list()
            Muscarinic_true = list()

            for seci in new_sim.neuron.h.allsec():

                for segi in seci:

                    if str('Unknown[1]') in str(segi) and str('axon') not in str(segi) and Muscarinic==True:

                        print('Muscarinic intracellular')
                        print(segi)
                        ach_conce = new_sim.neuron.h.concACh(segi)
                        pointer_concentration=ach_conce._ref_concentration
                        Minput = new_sim.neuron.h.M4(segi)
                        new_sim.neuron.h.setpointer(pointer_concentration,'conc_ACH',Minput)
                        pointer_modulation=Minput._ref_Ach_M4R
                        new_sim.neuron.h.setpointer(pointer_modulation,'muscarinic_modulation', segi.kir23_lts)
                        Muscarinic_true.append([ach_conce,Minput,pointer_modulation,pointer_concentration])
                        muscarinic_recording= new_sim.neuron.h.Vector()

                        if "dend" in str(segi):
                                segi.correction_kir23_lts = 1
                                #import pdb
                                #pdb.set_trace()
                        if "soma" in str(segi):
                                segi.correction_kir23_lts = 1
                                #import pdb
                                #pdb.set_trace()
                        muscarinic_recording.record(Minput._ref_Ach_M4R)
                        acetyl_recording= new_sim.neuron.h.Vector()
                        acetyl_recording.record(ach_conce._ref_concentration)

            if Muscarinic==True:
                for syn in Muscarinic_true:
                    print("Adding muscarinic")
                    nc_musc = new_sim.neuron.h.NetCon(getattr(seg,'_ref_v'),syn[0])
                    nc_musc.delay=0
                    nc_musc.threshold=0
                    nc_musc.weight[0]=parameters_file["muscarinic_weight"]
                    muscsynapses.append([nc_musc,syn[0]])


        elif str(Cells['ChIN_1']['soma_access']) in str(seg):

            #print(str(Cells['ChIN_1']['soma_access']))
            print('ATTACHING LTS SPIKING TO THE ChIN')

            LTSHzChecker = new_sim.neuron.h.IntFire1()
            LTSHzChecker.tau=1000
            LTSHzChecker.refrac = 10
            Int_LTS=new_sim.neuron.h.Vector()
            Int_LTS.record(LTSHzChecker._ref_m)

            nc_LTS = neuron.h.NetCon(getattr(seg,'_ref_v'),LTSHzChecker)
            nc_LTS.weight[0] = Checkerweight[1]
            nc_LTS.delay =0
            nc_LTS.threshold =-40


            Nitric_oxide = list()

            for seci in new_sim.neuron.h.allsec():

                for segi in seci:

                    if NO==True and str(Cells['ChIN_0']['soma_access']) in str(segi):
                        print("NO")          
                        NOinput = new_sim.neuron.h.NO(segi)
                        Nitric_oxide.append(NOinput)


            if NO==True:
                for syn in Nitric_oxide:
                    print("Nitric oxide")
                    nc_no = new_sim.neuron.h.NetCon(LTSHzChecker,syn)
                    nc_no.delay=0
                    nc_no.threshold=0
                    nc_no.weight[0]=parameters_file["nitric_oxide_weight"]
                    NOsynapses.append(nc_no)


for sec in new_sim.neuron.h.allsec():
    if "" in sec.name():
        print(new_sim.neuron.h.psection())

tSave = new_sim.neuron.h.Vector()
tSave.record(new_sim.neuron.h._ref_t)

print("Input to network")
print("----------------")
print(random_stream_Crtx_LTS)
print("---------------- \n")

print(Synapse_TOTAL)


print(random_stream_TH)
print("Simulating whole network \n")
print(str(configs) + '\n')


new_sim.neuron.h.tstop=parameters_file["simulation_time"]
new_sim.neuron.h.run()


k=0

name_cells={"ChIN_0":"ChIN","ChIN_1":"LTS"}

for cell,cell_info in Cells.items():
    print(cell)
    cell_v=cell_info["soma_voltage"]    
    plt.figure(k)
    plt.plot(tSave,cell_v,label=name_cells[cell],c='black')
    plt.title("ChIN_LTS")
    plt.legend()
    #plt.show() 
    plt.savefig('Output/Fig-7/cortical/'+parameters_file["single-cell-models"][0].split("/")[-1].split(".")[0]+'-cortex'+str(cortex_frequency)+'-thalamus-'+str(thalamus_frequency)+'-checker-'+str(checkerweight)+str(Synapse_TOTAL)+str(name_cells[cell])+'.svg')
    np.savetxt('Output/Fig-7/cortical/Opt0-voltage-cortex-'+str(cortex_frequency)+'-thalamus-'+str(thalamus_frequency)+'-checker-'+str(checkerweight)+str(Synapse_TOTAL)+name_cells[cell]+'.txt',[tSave,cell_v])
    plt.clf()

    k=k+1
  
plt.figure(2)
plt.plot(tSave,Int_LTS,label='LTS')
plt.title('Checker on LTS')
plt.legend()
plt.savefig('Output/Fig-7/cortical/'+parameters_file["single-cell-models"][0].split("/")[-1].split(".")[0]+'-Checker-LTS-cortex-'+str(cortex_frequency)+'-thalamus-'+str(thalamus_frequency)+'-checker-'+str(checkerweight)+str(Synapse_TOTAL)+'.svg')
plt.clf()

plt.figure(4)
plt.plot(tSave,muscarinic_recording)
plt.title('Muscarinic recording')
plt.savefig('Output/Fig-7/cortical/'+parameters_file["single-cell-models"][0].split("/")[-1].split(".")[0]+'-M4-recording-cortex-'+str(cortex_frequency)+'-thalamus-'+str(thalamus_frequency)+'-checker-'+str(checkerweight)+str(Synapse_TOTAL)+'.svg')
plt.clf()

plt.figure(5)
plt.plot(tSave,acetyl_recording)
plt.title('Acetylcholine activation')
plt.savefig('Output/Fig-7/cortical/'+parameters_file["single-cell-models"][0].split("/")[-1].split(".")[0]+'-Acetylcholine-cortex-'+str(cortex_frequency)+'-thalamus-'+str(thalamus_frequency)+'-checker-'+str(checkerweight)+str(Synapse_TOTAL)+'.svg')
plt.clf()



print('wall time: {}s'.format(time.time() - start_time))



