"""
Multicompartment in the mouse brain.
This is the script for the 7 compartments Inulin test case
This scriipt cqn qlso be used to simulate the 4-compartment case (just put to zero the fluid transfer between blood and PVS compartments)


The equations are
1) equations for the pore pressures
2) advection-diffusion equations for the tracer concentration within the pores


For this script we consider 4 compartments:
 - interstitial space (0)
 - arteries (1)
 - veins (2)
 - capillaries (3)
 - PVS arteries (4)
 - PVS veins (5)
 - PVS capillaries (6)

Units are: mm for space and second for time
Standard use for this script is

python3 name_script.py

"""
from dolfin import *
import numpy as np
from ts_storage import TimeSeriesStorage
from pathlib import Path
import sys
import argparse
from typing import Sequence, Optional


def ForwardModel(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-ecs_factor", type=float,
                        help="Scaling factor for ECS")
    parser.add_argument("-pvs_cap_factor", type=float,
                        help="Scaling factor for PVS_CAP")
    parser.add_argument("-pCFS_factor", type=float,
                        help="Scaling factor for pCFS")
    parser.add_argument("-pvsc_ECS_transfer_factor",
                        type=float, help="Scaling factor for pvsc_ECS")

    args = vars(parser.parse_args(argv))
    forward(**args)
    return 0


def forward(ecs_factor: float, pvs_cap_factor: float, pCFS_factor: float, pvsc_ECS_transfer_factor: float):

    # Choose Boundary type
    #BC_type = "Homogeneous"
    # BC_type = "Conservation"  # You have three choices: Conservation, Homogeneous, Decay
    BC_type = "Decay"

    # Temporal parameters
    dt = 200.
    T = 3600.*6
    decay = 0.01/60.

    # Load mesh
    res = 16
    meshfile = Path(f"mesh/mesh{res}.h5")
    mesh = Mesh()
    hdf = HDF5File(mesh.mpi_comm(), str(meshfile), "r")
    hdf.read(mesh, "/mesh", False)
    SD = MeshFunction("size_t", mesh, mesh.topology().dim())
    hdf.read(SD, "/subdomains")
    bnd = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    hdf.read(bnd, "/boundaries")

    # define compartments
    comp = ['ISF', 'arteries', 'veins', 'capillaries',
            'PVS arteries', 'PVS veins', 'PVS capillaries']
    ncomp = len(comp)
    geo = mesh.ufl_cell()
    h = mesh.hmin()
    # print(h)

    # Finite element functions (P1 for everything)
    P1 = FiniteElement('CG', geo, 1)
    P2 = FiniteElement('Lagrange', geo, 2)
    ME = MixedElement(ncomp*[P2])
    Q = FunctionSpace(mesh, ME)
    V = VectorFunctionSpace(mesh, 'Lagrange', 2)
    VV = FunctionSpace(mesh, P2)
    # print(Q.dim())
    p = TrialFunctions(Q)
    q = TestFunctions(Q)

    # Load measure data
    dx = Measure("dx", domain=mesh, subdomain_data=SD)  # Volume
    ds = Measure('ds')()  # surface

    # Compute surface area of the mesh and its volume
    surface_area = assemble(1.0*ds(domain=mesh))  # in mm^2
    #print('Surface area: ', surface_area, ' mm^2')
    brain_volume = assemble(1*dx(domain=mesh))
    #print('brain volume: ', brain_volume, ' mm^3')
    Vcsf = 0.1 * brain_volume
    #print('brain volume: ', Vcsf, ' mm^3')
    n = FacetNormal(mesh)

    # PHYSICAL PARAMETERS
    # porosity
    phi0 = ncomp*[Constant(5e-8)]  # Porosity V_i/V_Total
    phi0[0] = 0.14
    phi0[1] = 0.00658
    phi0[2] = 0.02303
    phi0[3] = 0.00329
    phi0[4] = 6.e-4
    phi0[5] = 2.1e-3
    phi0[6] = 3.e-4

    # viscosity
    nu = np.array([0.0]*ncomp)
    nu[0] = 7.0e-4
    nu[1] = 2.67e-3  # blood viscosity
    nu[2] = 2.67e-3
    nu[3] = 2.67e-3
    nu[4] = 7.0e-4
    nu[5] = 7.0e-4
    nu[6] = 7.0e-4

    # Permeability of fluid
    kappa_f = np.array([0.0]*ncomp)
    kappa_f[0] = 2.0e-11*ecs_factor
    kappa_f[1] = 3.29478e-06
    kappa_f[2] = 6.58956e-06
    kappa_f[3] = 1.14276e-09
    kappa_f[4] = 1.0e-11*pvs_cap_factor
    kappa_f[5] = 6.514285714285714e-09
    kappa_f[6] = 3.5359801488833745e-13

    # transfer coefficients

    w_apa = 5.89e-09
    w_vpv = 1.26e-10
    w_cpc = 2.98e-09

    # WARNING: Comment the floowing 3 lines to have 7 compartments
    # w_apa = 0
    # w_vpv = 0
    # w_cpc = 0

    w_pae = 2.1932017763179826e-07  # *pvsc_ECS_transfer_factor
    w_pve = 1.9533203320332029e-07
    w_pce = 9.976258977745193e-10

    w_ac = 1.05010501050105e-07
    w_cv = 5.250525052505252e-07
    w_papc = 2.5002500250025003e-08
    w_pcpv = 1.0001000100010001e-07

    # Fluid pressure exchange
    gamma = np.array([[0, 0, 0, 0, w_pae, w_pve, w_pce],
                      [0, 0, 0, w_ac, w_apa, 0, 0],
                      [0, 0, 0, w_cv, 0, w_vpv, 0],
                      [0, w_ac, w_cv, 0, 0, 0, w_cpc],
                      [w_pae, w_apa, 0, 0, 0, 0, w_papc],
                      [w_pve, 0, w_vpv, 0, 0, 0, w_pcpv],
                      [w_pce, 0, 0, w_cpc, w_papc, w_pcpv, 0]])
    # print(gamma)

    osmo_cap = 20.0*133.33
    # osmo_e = 0.*osmo_cap  # If you consider the blood compartments, change this line to:
    osmo_e = 0.2*osmo_cap
    osmo = np.array([osmo_e, osmo_cap, osmo_cap,
                    osmo_cap, osmo_e, osmo_e, osmo_e])

    # INITIAL CONDITIONS

    F = 0
    # Solving steady state pressure
    for i in range(ncomp):
        F += kappa_f[i]/nu[i]*inner(grad(p[i]), grad(q[i]))*dx

        for j in range(ncomp):
            if i != j:
                F -= Constant(gamma[i][j])*inner(p[j]-p[i] -
                                                 (osmo[j]-osmo[i]), q[i])*dx

    # Boundary conditions
    F -= 1/(0.01*133.32*60.0/1000.0*1.0e6) * \
        (120.0*133.33*q[1]*ds - p[1]*q[1]*ds)
    F -= (3.13e-7)*(3.26*133.33*q[0]*ds - p[0]*q[0]*ds)

    p_CSF = 4.74*133.333*pCSF_factor
    F -= (1.25e-6)*(p_CSF*q[4]*ds - p[4]*q[4]*ds)

    bcs = [DirichletBC(Q.sub(2), Constant(7.0*133.33), 'on_boundary'),
           DirichletBC(Q.sub(5), Constant(3.26*133.33), 'on_boundary')]

    # Solving pressure equations
    print("Solving pressure system...", sep="")
    A = assemble(lhs(F))
    b = assemble(rhs(F))
    [bc.apply(A, b) for bc in bcs]
    p_ = Function(Q)
    solve(A, p_.vector(), b, 'gmres', 'ilu')
    print("Done.")

    # Some manips to get the pressure fields for the pressure equations
    p_new = p_.split(True)
    p1, p2, p3, p4, p5, p6, p7 = p_.split(True)
    p_0 = Function(VV)
    assign(p_0, p_new[0])
    p_1 = Function(VV)
    assign(p_1, p_new[4])
    p_2 = Function(VV)
    assign(p_2, p_new[5])
    p_3 = Function(VV)
    assign(p_3, p_new[6])

    velocities = np.array([[float(1/brain_volume*assemble(kappa_f[i]/nu[i] *
                                                          dot(grad(p_new[i]), grad(p_new[i]))**0.5*dx)) for i in [0, 4, 5, 6]], ])

    ########
    # PART 2: the diffusion-convection equations
    ########

    ecs_permeability = float(kappa_f[6])
    # Go to 4 compartments because we do not have transport in blood
    comp = ['IS', 'PVS arteries', 'PVS veins', 'PVS capillaries']
    ncomp = len(comp)

    # Diffusion Coefficient
    D_free = 2.98e-4
    D_eff = 1.03e-4
    # porosity

    phi0 = np.append(phi0[0], phi0[4:])

    # viscosity
    nu = np.append(nu[0], nu[4:])
    # Permeability of fluid
    kappa_f = np.append(kappa_f[0], kappa_f[4:])

    # INULIN exchange
    sigma_reflect_AEF = 0.2
    g_ae = w_pae*(1-sigma_reflect_AEF)
    g_ce = w_pce*(1-sigma_reflect_AEF)
    g_ve = w_pve*(1-sigma_reflect_AEF)
    g_ac = w_papc
    g_cv = w_pcpv
    gamma_tilde = np.array([[0, g_ae, g_ve, g_ce], [g_ae, 0, 0, g_ac], [
        g_ve, 0, 0, g_cv], [g_ce, g_ac, g_cv, 0]])

    # From diffusion
    l_ae = 3.744331208456129e-05*pvsc_ECS_transfer_factor
    l_ce = 1.7198523872626686e-05
    l_ve = 3.744331208456129e-05
    l_ac = 0
    l_cv = 0
    lmbd = np.array([[0, l_ae, l_ve, l_ce], [l_ae, 0, 0, l_ac],
                    [l_ve, 0, 0, l_cv], [l_ce, l_ac, l_cv, 0]])

    ME = MixedElement(ncomp*[P2])
    Q = FunctionSpace(mesh, ME)
    p_new = Function(Q)
    assign(p_new, [p_0, p_1, p_2, p_3])

    # CONCENTRATIONS EQUATIONS
    P1 = FiniteElement('CG', geo, 1)
    P2 = FiniteElement('Lagrange', geo, 2)
    ME = MixedElement(ncomp*[P1])
    Q = FunctionSpace(mesh, ME)
    V = VectorFunctionSpace(mesh, 'Lagrange', 1)
    VV = FunctionSpace(mesh, P1)
    # print(Q.dim())
    p = TrialFunctions(Q)
    q = TestFunctions(Q)
    # Gaussian initial condition.
    #print("Computing initial conditions")
    center = (4., 2., 3.)
    spread = 1.0
    u0 = Expression(
        "exp(- (pow(x[0]-s[0], 2) + pow(x[1]-s[1], 2) + pow(x[2]-s[2], 2)) / (b * b))",
        degree=1, b=Constant(spread), s=Constant(center)
    )

    # Project u0 to have a homogeneous Dirichlet boundary (might exists a lot better approaches)
    V = FunctionSpace(mesh, "Lagrange", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    a0 = u * v * dx
    L0 = u0 * v * dx
    init_c_ecs = Function(V)
    solve(a0 == L0, init_c_ecs, bcs=[
          DirichletBC(V, Constant(0.), "on_boundary")])
    init_c_ecs = project(0.14*init_c_ecs/(phi0[0]), V)

    G = 0  # init variational formsplit
    # function to save concentrations at the previous time step
    cn = Function(Q)

    # Variational form for the tracer concentration equations

    def Max(a, b): return (a+b+abs(a-b))/Constant(2)
    def Min(a, b): return (a+b-abs(a-b))/Constant(2)

    for i in range(ncomp):
        if i == 0:  # If in ISF grey and white matter diffusion
            G += Constant(D_eff)*inner(grad(p[i]), grad(q[i]))*dx
        else:  # Otherwise free diffusion
            G += Constant(D_eff)*inner(grad(p[i]), grad(q[i]))*dx

        G += 1/dt*inner(p[i]-cn[i], q[i])*dx

        G += Constant(kappa_f[i]/(phi0[i]*nu[i])) * \
            inner(p[i], inner(grad(p_new[i]), grad(q[i])))*dx
        # mass exchange
        for j in range(ncomp):
            if i != j:
                G += Constant(lmbd[i][j])*inner(p[i]-p[j], q[i])*1./phi0[i]*dx
                #G += Constant(gamma_tilde[i][j])*Max(p_new[i]-p_new[j],0)*p[i]*q[i]*1./phi0[i]*dx
                #G += Constant(gamma_tilde[i][j])*Min(p_new[i]-p_new[j],0)*p[j]*q[i]*1./phi0[i]*dx
                G += Constant(gamma_tilde[i][j])*(p_new[i] -
                                                  p_new[j])*(p[j]+p[i])*0.5*q[i]*1./phi0[i]*dx

    # Boundary conditions for tracer concentration
    if BC_type == "Homogeneous":
        A_in = Expression('C', C=0, degree=1)
        bcs_conc = [DirichletBC(Q.sub(1), A_in, 'on_boundary'), DirichletBC(
            Q.sub(2), A_in, 'on_boundary'), DirichletBC(Q.sub(0), A_in, 'on_boundary')]
    elif BC_type == "Conservation":
        g0 = 0.0  # Zero concentration at the beginning
        g = Constant(g0)
        g1 = Constant(g)
        g2 = Constant(g)
        g00 = Constant(g)
        bcs_conc = [DirichletBC(Q.sub(1), g1, 'on_boundary'), DirichletBC(
            Q.sub(2), g2, 'on_boundary'), DirichletBC(Q.sub(0), g00, 'on_boundary')]
    elif BC_type == "Decay":
        g0 = 0.0  # Zero concentration at the beginning
        g = Constant(g0)
        g1 = Constant(g)
        g2 = Constant(g)
        g00 = Constant(g)
        bcs_conc = [DirichletBC(Q.sub(1), g1, 'on_boundary'), DirichletBC(
            Q.sub(2), g2, 'on_boundary'), DirichletBC(Q.sub(0), g00, 'on_boundary')]
    elif BC_type == "zeroNeum":
        bcs_conc = []
    else:
        print('Wrong Boundary conditions')
        exit(0)

    # ASSEMBLING
    #print("Assembling diffusion problem")
    a_conc = lhs(G)
    L_conc = rhs(G)
    A_conc = assemble(a_conc)

    czero = Expression('0.0', degree=1)
    #init_c_ecs = interpolate(u0, Q.sub(0).collapse())
    init_other = interpolate(czero, V)
    c_ = Function(Q)  # Function to save solution
    assign(c_, [init_c_ecs, init_other, init_other, init_other])
    [bc.apply(c_.vector()) for bc in bcs_conc]
    c1 = c_.split(True)

    results_path = Path(
        f"results_sensitivity/results-{BC_type}-mesh{res}-dt{int(dt)}-inulin-{len(comp)}comps-kecs{ecs_factor:.4f}-kpvscap{pvs_cap_factor:.4f}-pCSF{pCSF_factor:.4f}-wPCE{pvsc_ECS_transfer_factor:.4f}")
    results_path.mkdir(parents=True, exist_ok=True)

    storage_cecs = TimeSeriesStorage(
        "w", results_path, mesh=mesh, V=VV, name="ecs")
    storage_carteries = TimeSeriesStorage(
        "w", results_path, mesh=mesh, V=VV, name="arteries")
    storage_cveins = TimeSeriesStorage(
        "w", results_path, mesh=mesh, V=VV, name="veins")
    storage_ccap = TimeSeriesStorage(
        "w", results_path, mesh=mesh, V=VV, name="capillaries")

    init_c_ecs.rename('concentration', '')
    storage_cecs.write(c1[0], 0.)
    storage_carteries.write(c1[1], 0.)
    storage_cveins.write(c1[2], 0.)
    storage_ccap.write(c1[3], 0.)

    # Store some vector values for total inulin amount in brain
    N0 = Constant(phi0[0] * assemble(init_c_ecs * dx))
    amount = np.zeros((int(T / dt + 1), ncomp+1))
    amount[0, 1] = 1.0

    # Storage for point concentration
    concentration_p = np.zeros_like(amount)
    concentration_p[0] = init_c_ecs(center)

    # Total tracer amount in system
    N = np.zeros_like(amount)
    N[0] = N0
    #print("mass_init =" +str(N[0]))
    times = np.zeros(len(amount))
    t = 0.0
    t += dt
    it = 1
    transfer_arteries = np.zeros_like(amount)
    print("solving diffusion equations")
    # Time steping
    while t < T + dt/2:  # Added dt/2 to ensure final time included.
        #print('t = ', t)
        cn.assign(c_)

        b_conc = assemble(L_conc)
        [bc.apply(A_conc, b_conc) for bc in bcs_conc]

        # Solve
        #print("solving diffusion equations")

        solve(A_conc, c_.vector(), b_conc, 'gmres', 'ilu')
        c1 = c_.split(True)
        storage_cecs.write(c1[0], t)  # store the ISF concentration
        storage_carteries.write(c1[1], t)  # store the arterial concentration
        storage_ccap.write(c1[3], t)  # store the capillary concentration
        mass_in = 0
        for j in range(ncomp):
            mass_in += assemble(phi0[j]*c1[j]*dx)
        # storage_ccap.write(c1[3], t) # store the capillary concentration
        ce_ = assemble(phi0[0]*c1[0]*dx)/float(N0)
        ca_ = assemble(phi0[1]*c1[1]*dx)/float(N0)
        cv_ = assemble(phi0[2]*c1[2]*dx)/float(N0)
        cc_ = assemble(phi0[3]*c1[3]*dx)/float(N0)
        amount[it, 1:] = np.array([ce_, ca_, cv_, cc_])
        amount[it, 0] = t/60
        times[it] = t/60
        concentration_p[it] = c1[0](center)

        # Prepare boundary conditions
        if BC_type == "Conservation":
            mass_out = 0
            mass_in = 0
            mmass = 0
            for j in range(ncomp):
                mass_in += assemble(phi0[j]*c_[j]*dx)
                mmass += assemble(phi0[j]*(c_[j] - cn[j])/dt * dx)
            transfer = 0.0
            for i in range(ncomp):
                for j in range(ncomp):
                    if i != j:
                        transfer += assemble(lmbd[i][j]*(c_[i]-c_[j])*dx)
                        transfer += assemble(gamma_tilde[i][j] *
                                             (p_new[i]-p_new[j])*(c_[j]+c_[i])*0.5*dx)

            #print("real mass total = " + str(amount[0]))
            #print("mass inside = " + str(mass_in))
            #print("g = " + str(float(g)))
            #print("Total transfer = " + str(transfer))
            g.assign(g - dt*mmass/Vcsf)

            g1.assign(g)
            g2.assign(g)
            g00.assign(g)

        elif BC_type == "Decay":
            #mass_out = assemble(phi0[0]*D_eff/Vcsf* grad(c1[0])*n * ds + phi0[0]*kappa_f[0]/(Vcsf*nu[0])*c1[0]*grad(p_new[0])*n *ds )
            mass_in = 0
            mmass = 0
            for j in range(ncomp):
                mass_in += assemble(phi0[j]*c_[j]*dx)
                mmass += assemble(phi0[j]*(c_[j] - cn[j])/dt * dx)
            g.assign((1/(1 + decay*dt))*(g - dt/Vcsf * (mmass)))
            g1.assign(g)
            g2.assign(g)
            g00.assign(g)

        it += 1

        t += dt

    # Storage
    storage_cecs.write(c1[0], t)  # store the ISF concentration
    storage_carteries.write(c1[1], t)  # store the arterial concentration
    storage_cveins.write(c1[2], t)  # store the venous concentration
    storage_ccap.write(c1[3], t)  # store the capillary concentration

    storage_cecs.store_info()
    storage_carteries.store_info()
    storage_cveins.store_info()
    storage_ccap.store_info()

    storage_cecs.close()
    storage_carteries.close()
    storage_cveins.close()
    storage_ccap.close()

    """
    # Compare clearence with diffusion
    plt.plot(times, amount, label="Multicompartment")
    plt.legend()
    # plt.ylim(0, None, auto=True)
    plt.savefig(results_path / f"inulin-multicomp-{BC_type}-dt{dt}-res{res}-{len(comp)}comps.png", bbox_inches='tight')
    plt.show()
    # save the clearance
    """
    np.savetxt(results_path / f"amount-multicomp-{BC_type}-dt{dt}-res{res}-{len(comp)}comps.csv",
               amount, delimiter=",", header='t,ecs,pa,pv,pc', comments='')
    np.savetxt(results_path / f"velocities-multicomp-{BC_type}-dt{dt}-res{res}-{len(comp)}comps.csv",
               velocities, delimiter=",", header='ecs,pa,pv,pc', comments='')


if __name__ == "__main__":
    ForwardModel()
