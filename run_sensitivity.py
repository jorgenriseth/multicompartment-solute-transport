import argparse
from pylab import *
from sensitivity import ForwardModel
from typing import List
import pandas as pd


def run_sensitivity(ecs_factor: List[float], pvs_cap_factor: List[float],
                    pCSF_factor: List[float], pvsc_ECS_factor: List[float]):
    sim_count = 1
    N = len(ecs_factor)*len(pvs_cap_factor)*len(pCSF_factor)*len(pvsc_ECS_factor)

    for k_ecs in ecs_factor:
        for k_pvs in pvs_cap_factor:
            for p in pCSF_factor:
                for g in pvsc_ECS_factor:
                    print(f"Running simulation {sim_count} of of {N}")
                    args = ["-ecs_factor", str(k_ecs), "-pvs_cap_factor", str(k_pvs),
                            "-pCSF_factor", str(p), "-pvsc_ECS_transfer_factor", str(g)]
                    ForwardModel(args)
                    sim_count += 1


def plot_sensitivity(ecs_factor: List[float], pvs_cap_factor: List[float],
                     pCSF_factor: List[float], pvsc_ECS_factor: List[float]):
    ecs = []*len(ecs_factor)
    subplot(2, 2, 1)
    k_pvs = pvs_cap_factor[0]
    p = pCSF_factor[0]
    g = pvsc_ECS_factor[0]

    for k in ecs_factor:
        results = f'results_sensitivity/results-Decay-mesh16-dt200-inulin-4comps-kecs{k:.4f}-kpvscap{k_pvs:.4f}-pCSF{p:.4f}-wPCE{g:.4f}/amount-multicomp-Decay-dt200.0-res16-4comps.csv'
        velocity = f'results_sensitivity/results-Decay-mesh16-dt200-inulin-4comps-kecs{k:.4f}-kpvscap{k_pvs:.4f}-pCSF{p:.4f}-wPCE{g:.4f}/velocities-multicomp-Decay-dt200.0-res16-4comps.csv'
        A = pd.read_csv(results)
        V = pd.read_csv(velocity)
        print('vel: ECS = %.2e, pvs_a = %.2e, pvs_c = %.2e, pvs_v = %.2e' %
              (V.ecs[0], V.pa[0], V.pc[0], V.pv[0]))
        plot(A.ecs)

    xlabel('Time [min]')
    ylabel('Relative concentration')
    legend(ecs_factor)
    title('ECS permeability')
    subplot(2, 2, 2)
    for k_pvs in pvs_cap_factor:
        results = f'results_sensitivity/results-Decay-mesh16-dt200-inulin-4comps-kecs{k:.4f}-kpvscap{k_pvs:.4f}-pCSF{p:.4f}-wPCE{g:.4f}/amount-multicomp-Decay-dt200.0-res16-4comps.csv'
        velocity = f'results_sensitivity/results-Decay-mesh16-dt200-inulin-4comps-kecs{k:.4f}-kpvscap{k_pvs:.4f}-pCSF{p:.4f}-wPCE{g:.4f}/velocities-multicomp-Decay-dt200.0-res16-4comps.csv'
        A = pd.read_csv(results)
        V = pd.read_csv(velocity)
        print('vel: ECS = %.2e, pvs_a = %.2e, pvs_c = %.2e, pvs_v = %.2e' %
              (V.ecs[0], V.pa[0], V.pc[0], V.pv[0]))
        plot(A.ecs)
    xlabel('Time [min]')
    ylabel('Relative concentration')
    legend(pvs_cap_factor)
    title(r'PVS$_a$ permeability')
    subplot(2, 2, 3)
    for p in pCSF_factor:
        results = f'results_sensitivity/results-Decay-mesh16-dt200-inulin-4comps-kecs{k:.4f}-kpvscap{k_pvs:.4f}-pCSF{p:.4f}-wPCE{g:.4f}/amount-multicomp-Decay-dt200.0-res16-4comps.csv'
        velocity = f'results_sensitivity/results-Decay-mesh16-dt200-inulin-4comps-kecs{k:.4f}-kpvscap{k_pvs:.4f}-pCSF{p:.4f}-wPCE{g:.4f}/velocities-multicomp-Decay-dt200.0-res16-4comps.csv'
        A = pd.read_csv(results)
        V = pd.read_csv(velocity)
        print('vel: ECS = %.2e, pvs_a = %.2e, pvs_c = %.2e, pvs_v = %.2e' %
              (V.ecs[0], V.pa[0], V.pc[0], V.pv[0]))
        plot(A.ecs, label=f'CSFp - {p}')
    xlabel('Time [min]')
    ylabel('Relative concentration')
    title('CSF pressure')
    legend(pCSF_factor)
    subplot(2, 2, 4)
    for g in pvsc_ECS_factor:
        results = f'results_sensitivity/results-Decay-mesh16-dt200-inulin-4comps-kecs{k:.4f}-kpvscap{k_pvs:.4f}-pCSF{p:.4f}-wPCE{g:.4f}/amount-multicomp-Decay-dt200.0-res16-4comps.csv'
        velocity = f'results_sensitivity/results-Decay-mesh16-dt200-inulin-4comps-kecs{k:.4f}-kpvscap{k_pvs:.4f}-pCSF{p:.4f}-wPCE{g:.4f}/velocities-multicomp-Decay-dt200.0-res16-4comps.csv'
        A = pd.read_csv(results)
        V = pd.read_csv(velocity)
        print('vel: ECS = %.2e, pvs_a = %.2e, pvs_c = %.2e, pvs_v = %.2e' %
              (V.ecs[0], V.pa[0], V.pc[0], V.pv[0]))
        plot(A.ecs, label=f'pvs-ecs - {g}')
    xlabel('Time [min]')
    ylabel('Relative concentration')
    title(r'PVS$_a$-ECS transfer')
    legend(pvsc_ECS_factor)
    tight_layout()
    figure()

    subplot(2, 2, 1)

    show()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-cmd", type=str, choices=["plot", "run"],
                        help="Run sensitivity model or plot results", required=True)
    args = vars(parser.parse_args())

    ecs_factor = [0.1, 1.0, 10, 100, 1000]
    # In the last run the runs with 1,1,1,pvs_cap_factor represents pvs_a_permeability factor. For combinations = still pcs_cap factor.
    pvs_cap_factor = [0.01, 0.1, 1.0, 10, 100]
    pCSF_factor = [0.5, 0.75, 1.0, 1.5, 2]

    # THIS IS NOW DIFFUSIVE EXCHANGE BETWEEN PERIARTERIAL AND ECS
    pvsc_ECS_factor = [0.01, 0.1, 1.0, 10, 100]

    ecs_factor, pvs_cap_factor, pvsc_ECS_factor = [1.0], [1.0], [1.0]

    if args["cmd"] == "plot":
        plot_sensitivity()
    elif args["cmd"] == "run":
        run_sensitivity(ecs_factor, pvs_cap_factor, pCSF_factor, pvsc_ECS_factor)
    else:
        ValueError("Invalid command")
