import os
import sys
from IPython import embed
from pylab import *

k_e_factor = [0.1, 1.0, 10, 100, 1000]
k_pa_factor = [0.01, 0.1, 1, 10, 100]
k_pc_factor = [0.01, 0.1, 1.0, 10, 100]
pCSF_factor = [0.5, 0.75, 1.0, 1.5, 2]
p_ae_diff_factor = [0.01, 0.1, 1.0, 10, 100]
gamma_paa_factor = [0.01, 0.1, 1.0, 10, 100]
D_factor = [0.5, 0.75, 1.0, 1.5, 2]
phi_pa_factor = [0.1, 1, 10, 20, 50]

#k_e_factor, p_pa_factor, k_pc_factor, pCSF_factor, p_ae_diff_factor, gamma_paa_factor, D_factor = 7*[[1]]

#ecs_factor, pvs_cap_factor, pvsc_ECS_factor = [1.0],[1.0],[1.0]


def run_sensitivity(kind='serial'):
    sim_count = 1
    
    
    if kind != 'serial':
        params = [k_e_factor, k_pa_factor, k_pc_factor, pCSF_factor, p_ae_diff_factor, gamma_paa_factor, D_factor, phi_pa_factor]
        N = sum([len(i) for i in params])
        for k_ecs in k_e_factor:
            for k_pa in k_pa_factor:
                for k_pvc in k_pc_factor:
                    for p in pCSF_factor:
                        for pae_factor in p_ae_diff_factor:                                 
                            for g in pvsc_ECS_factor:
                                print(f"Running simulation {sim_count} of of {N}")
                                os.system(f'python3 sensitivity.py {k_ecs} {k_pvs} {p} {g}')
                                sim_count += 1
    else:
        params = [k_e_factor, k_pa_factor, k_pc_factor, pCSF_factor, p_ae_diff_factor, gamma_paa_factor, D_factor, phi_pa_factor]
        N = sum([len(i) for i in params])
        n = len(params)
        i = 0
        sim_count = 1
        for factors in params:
            i += 1
            if factors != [1]:
                for p in factors:
                    print(f"Running simulation {sim_count} of of {N}")
                    cmd_args = n*'1 '
                    cmd_args = cmd_args[:2*(i-1)] + '%s '%p + cmd_args[2*i:]
                    cmd = f'python3 sensitivity.py ' + cmd_args
                    print(cmd)
                    os.system(cmd)
                    sim_count += 1


def plot_sensitivity(comps = 4):
    import pandas as pd

    params = [k_e_factor, k_pa_factor, k_pc_factor, pCSF_factor, p_ae_diff_factor, gamma_paa_factor, D_factor, phi_pa_factor]
    titles = [r'$\kappa_e$', r'$\kappa_{pa}$', r'$\kappa_{pc}$', r'$p_{SAS}$', r'$P_{pa,e}$', r'$\gamma_{pa,a}$', r'$D$', r'$\phi_{pa}$']
    figure(figsize=(16,10))
    for i in range(len(params)):
        d = {'k_e':1, 'k_pa':1, 'k_pc':1, 'pCSF':1, 'p_ae_diff':1, 'gamma_paa':1, 'D':1, 'phi_pa':1}
        factors = params[i]
        for p in factors:
            param = list(d.keys())[i]
            d[param] = p
            basefolder = 'results_sensitivity_%s_comp/results-Decay-mesh16-dt200-inulin-%scomps-'%(comps, comps)
            paramstring = '-'.join([f'{key}{value:.4f}' for key,value in d.items()])
            results = basefolder + paramstring + '/amount-multicomp-Decay-dt200.0-res16-4comps.csv'
            velocities = basefolder + paramstring + '/velocities-multicomp-Decay-dt200.0-res16-4comps.csv'
            #results = f"results_sensitivity/results-Decay-mesh16-dt200-inulin-4comps-k_e{d['k_e']:.4f}-k_pa{d['k_pa']:.4f}-k_pc{d['k_pc']:.4f}-pCSF{d['pCSF']:.4f}-p_ae_diff{d['p_ae_diff']:.4f}-gamma_paa{d['gamma_paa']:.4f}-D{d['D']:.4f}/amount-multicomp-Decay-dt200.0-res16-4comps.csv"
            #velocities = f"results_sensitivity/results-Decay-mesh16-dt200-inulin-4comps-k_e{d['k_e']:.4f}-k_pa{d['k_pa']:.4f}-k_pc{d['k_pc']:.4f}-pCSF{d['pCSF']:.4f}-p_ae_diff{d['p_ae_diff']:.4f}-gamma_paa{d['gamma_paa']:.4f}-D{d['D']:.4f}/velocities-multicomp-Decay-dt200.0-res16-4comps.csv"
            subplot(2,4,i+1)
            A = pd.read_csv(results)
            V = pd.read_csv(velocities)
            plot(A.ecs)
            print(titles[i], list(A['ecs'])[-1], float(V.ecs))
        legend(factors)
        title(titles[i])
             

    show()


cmd = sys.argv[1]
if cmd == 'plot':
    try:
        comps = int(sys.argv[2])
        plot_sensitivity(comps)
    except:    
        plot_sensitivity()
elif cmd == 'run':
    run_sensitivity() 

