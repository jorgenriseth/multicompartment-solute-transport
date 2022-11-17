import os
from IPython import embed
from pylab import *

ecs_factor = [0.1, 1.0, 10, 100, 1000]
pvs_cap_factor = [0.01, 0.1, 1.0, 10, 100]
pCSF_factor = [0.5, 0.75, 1.0, 1.5, 2]
pvsc_ECS_factor = [0.01, 0.1, 1.0, 10, 100]

#pvs_cap_factor, pCSF_factor, pvsc_ECS_factor = [1.0],[1.0],[1.0]


def download_sensitivity():
    sim_count = 1
    N = len(ecs_factor)*len(pvs_cap_factor)*len(pCSF_factor)*len(pvsc_ECS_factor)
    
    for k_ecs in ecs_factor:
        for k_pvs in pvs_cap_factor:
            for p in pCSF_factor:
                for g in pvsc_ECS_factor:
                    folder = f'results_sensitivity/results-Decay-mesh16-dt200-inulin-4comps-kecs{k_ecs:.4f}-kpvscap{k_pvs:.4f}-pCSF{p:.4f}-wPCE{g:.4f}'
                    results = folder + '/amount-multicomp-Decay-dt200.0-res16-4comps.csv'
                    os.system(f'mkdir {folder}')
                    os.system(f'scp -r vegard@calliope:Vegard/multicompartment-solute-transport/{folder}/*.csv {results}')
                    print(f'Downloading data set {sim_count} ot of {N}')
                    sim_count += 1
        
 
download_sensitivity()

