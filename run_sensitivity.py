import os
from IPython import embed
from pylab import *

ecs_factor = [0.1, 1.0, 10, 100, 1000]
pvs_cap_factor = [0.01, 0.1, 1.0, 10, 100]   # In the last run the runs with 1,1,1,pvs_cap_factor represents pvs_a_permeability factor. For combinations = still pcs_cap factor. 
pCSF_factor = [0.5, 0.75, 1.0, 1.5, 2]
pvsc_ECS_factor = [0.01, 0.1, 1.0, 10, 100]

#ecs_factor, pCSF_factor, pvsc_ECS_factor = [1.0],[1.0],[1.0]


def run_sensitivity():
    sim_count = 1
    N = len(ecs_factor)*len(pvs_cap_factor)*len(pCSF_factor)*len(pvsc_ECS_factor)
    
    for k_ecs in ecs_factor:
        for k_pvs in pvs_cap_factor:
            for p in pCSF_factor:
                for g in pvsc_ECS_factor:
                    print(f"Running simulation {sim_count} of of {N}")
                    os.system(f'python3 sensitivity.py {k_ecs} {k_pvs} {p} {g}')
                    sim_count += 1
        
 
k_pvs, p, g = 1.0,1.0,1.0
def plot_sensitivity():
    import pandas as pd
    ecs = []*len(ecs_factor)
    subplot(2,2,1)
    k_pvs, p, g = 1.0,1.0,1.0
    for k in ecs_factor:  
        results = f'results_sensitivity/results-Decay-mesh16-dt200-inulin-4comps-kecs{k:.4f}-kpvscap{k_pvs:.4f}-pCSF{p:.4f}-wPCE{g:.4f}/amount-multicomp-Decay-dt200.0-res16-4comps.csv'
        A = pd.read_csv(results)
        plot(A.ecs)
    legend(ecs_factor)
    title('ECS permeability')
    subplot(2,2,2)
    k, p, g = 1.0,1.0,1.0
    for k_pvs in pvs_cap_factor:  
        results = f'results_sensitivity/results-Decay-mesh16-dt200-inulin-4comps-kecs{k:.4f}-kpvscap{k_pvs:.4f}-pCSF{p:.4f}-wPCE{g:.4f}/amount-multicomp-Decay-dt200.0-res16-4comps.csv'
        A = pd.read_csv(results)
        plot(A.ecs)
    legend(pvs_cap_factor)
    title('PVS_a permeability')
    subplot(2,2,3)
    k, k_pvs, g = 1.0,1.0,1.0
    for p in pCSF_factor:  
        results = f'results_sensitivity/results-Decay-mesh16-dt200-inulin-4comps-kecs{k:.4f}-kpvscap{k_pvs:.4f}-pCSF{p:.4f}-wPCE{g:.4f}/amount-multicomp-Decay-dt200.0-res16-4comps.csv'
        A = pd.read_csv(results)
        plot(A.ecs,label=f'CSFp - {p}')
    title('CSF pressure')
    legend(pCSF_factor) 
    subplot(2,2,4)
    k, k_pvs, p = 1.0,1.0,1.0
    for g in pvsc_ECS_factor:  
        results = f'results_sensitivity/results-Decay-mesh16-dt200-inulin-4comps-kecs{k:.4f}-kpvscap{k_pvs:.4f}-pCSF{p:.4f}-wPCE{g:.4f}/amount-multicomp-Decay-dt200.0-res16-4comps.csv'
        A = pd.read_csv(results)
        plot(A.ecs,label=f'pvs-ecs - {g}')
    title('PVS_{cap}-ECS transfer')
    legend(pvsc_ECS_factor)
    show()
        
plot_sensitivity()

