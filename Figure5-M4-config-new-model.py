from Neuron_model_extended import NeuronModel
import sys
from neuron import h
import bluepyopt.ephys as ephys
import json     
import matplotlib.pyplot as plt
import numpy as np
import quantities as pq
import random
import time

with open("config-Figure5-M4.json",'r') as parameter_file:
        parameters_file = json.load(parameter_file)


new_sim=ephys.simulators.NrnSimulator(dt=0.025,cvode_active=False)

start_time=time.time()


Checkerweight=[0,0]

Cells = dict()

NOsynapses = list()

muscarinicsynapses=list()

thalamicsynapses = list()

corticalsynapses = list()

Gaussian = list()


configs=parameters_file["single-cell-models"]

name= ["ChIN","LTS"]

print()

for cells in range(2):

    with open(configs[cells],'r') as config_file:
        config = json.load(config_file)

    print(name[cells] + ' \n')

    morph = config["morphology"]
    print("morphology-DONE \n")
     
    param=config["parameters"]
    print("parameters-DONE \n")
        
    mech=config["mechanisms"]
    print("mechanisms DONE \n")

    info_cell = dict()

    model_cell=NeuronModel(param_file=param,morph_file=morph,mech_file=mech, cell_name=name[cells])

    new_sim=ephys.simulators.NrnSimulator(cvode_active=False)

    model_cell.instantiate(sim=new_sim)
    
    info_cell.update({"soma_sec": model_cell})

    info_cell.update({"soma_sim": new_sim})
    
    vSave = new_sim.neuron.h.Vector()
    
    for isec, sec in enumerate(model_cell.icell.soma):
        
        for seg in sec:
            if "0.5" in str(seg):
                gauss=new_sim.neuron.h.InGauss(seg)
                gauss.delay=0
                gauss.dur=10e10
                gauss.mean=0
                gauss.stdev=0.005
                Gaussian.append(gauss)
                vSave.record(getattr(seg,'_ref_v'))
                info_cell.update({"soma_access":seg})
                spike_time = new_sim.neuron.h.Vector()
                recording_netcon = new_sim.neuron.h.NetCon(getattr(seg,'_ref_v'),None, sec = sec)
                recording_netcon.threshold = 0
                recording_netcon.record(spike_time)
    
                                      
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
NO=False


for sec in new_sim.neuron.h.allsec():
    
    for seg in sec:
        
        if str(Cells['ChIN_0']['soma_access']) in str(seg):

            print('ATTACHING CHIN SPIKING TO THE LTS')

            ChINHzChecker = new_sim.neuron.h.IntFire1()
            ChINHzChecker.tau=1000
            ChINHzChecker.refrac = 10

            Int_ChIN=new_sim.neuron.h.Vector()
            Int_ChIN.record(ChINHzChecker._ref_m)
            
            nc_ChIN = neuron.h.NetCon(getattr(seg,'_ref_v'),ChINHzChecker)
            nc_ChIN.weight[0] = Checkerweight[0]
            nc_ChIN.delay = 0
            nc_ChIN.threshold = -40

            ChINprojectingsynapses = list()
            ChINnicotinicsynapses = list()
            Muscarinic_true = list()
            
            for seci in new_sim.neuron.h.allsec():
               
                for segi in seci:
   
                                    
                    if "LTS" in str(segi) and Muscarinic==True and "axon" not in str(segi):
                        #import pdb
                        #pdb.set_trace()                         
                        print('Muscarinic intracellular')

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
			
                
                            
            if Muscarinic==True:
                
                activity= parameters_file["muscarinic_activation_time"]

                VecStim_M=new_sim.neuron.h.VecStim()
                Vector_Activity=new_sim.neuron.h.Vector(activity)
                VecStim_M.play(Vector_Activity)
                muscarinicsynapses.append(VecStim_M)
                for syn in Muscarinic_true:
                    nc_musc = new_sim.neuron.h.NetCon(VecStim_M,syn[0])
                    nc_musc.delay=0
                    nc_musc.threshold=0
                    nc_musc.weight[0]=parameters_file["muscarinic_strength"]
                    muscarinicsynapses.append([nc_musc,syn[0]])
                        
            
        elif str(Cells['ChIN_1']['soma_access']) in str(seg):
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

                        NOinput = new_sim.neuron.h.NO(segi)
                        Nitric_oxide.append(NOinput)
                                              

            if NO==True:
                
                
                for syn in Nitric_oxide:
                    nc_no = new_sim.neuron.h.NetCon(LTSHzChecker,syn)
                    nc_no.delay=0
                    nc_no.threshold=0
                    nc_no.weight[0]=0
                    NOsynapses.append(nc_no)
                                  
        
                                
tSave = new_sim.neuron.h.Vector()
tSave.record(new_sim.neuron.h._ref_t)


print()
print("Simulating M4 input to LTS using \n")
print(str(configs) + '\n')

new_sim.neuron.h.tstop=parameters_file["simulation_time"]

new_sim.neuron.h.run()

print('wall time: {}s'.format(time.time() - start_time))

k=0


for cell,cell_info in Cells.items():
    
    print(cell)
    cell_v=cell_info["soma_voltage"]  

    if "ChIN_1" in cell:
            plt.figure(k)
            
            plt.plot(np.array(tSave),np.array(cell_v),label=cell,c='black')
            plt.savefig("Output/Fig-5-M4/Model_LTS-M4-"+parameters_file["single-cell-models"][1].split("/")[-1].split(".")[0]+'.svg')
            
            np.savetxt("Output/Fig-5-M4/Model_LTS-M4-"+parameters_file["single-cell-models"][1].split("/")[-1].split(".")[0]+'.txt',[tSave,cell_v])

    plt.title("ChIN_LTS")
    plt.legend()
    k=k+1

plt.show()



