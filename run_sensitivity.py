import os
import pandas as pd

ecs_factor = [0.1, 1, 10, 100, 1000]
pvs_cap_factor = [0.01, 0.1, 1, 10, 100]
pCSF_factor = [0.5, 0.75, 1, 1.5, 2]
pvsc_ECS_factor = [0.01, 0.1, 1, 10, 100]

k_pvs = 1
p = 1
g = 1


def run_sensitivity():
    for k_ecs in ecs_factor:
        os.system(f'python3 sensitivity.py {k_ecs} {k_pvs} {p} {g}')
        
 
def plot_sensitivity():
    ecs = []*len(ecs_factor)
    for k in ecs_factor:
        results = f'results_sensitivity/results-Decay-mesh16-dt600-inulin-4comps-kecs{k}-kpvscap1-pCSF1-wPCE1'
        pd.read_csv(results)
        
    embed()
        
run_sensitivity()
