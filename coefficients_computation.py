import numpy as np


def compute_coefficients(type_exp):
    """
        This function computes the coefficients for the multicompartment model.
        Input: none
        Output: dict of coefficients. The main script reads the coeffs from it.
    """
    # init main dict
    main_dict = {}

    # Small script to compute the right coefficients.
    volume_brain_human = 1.0e6 #mm^3
    volume_brain_rat = 2311.
    mmHgToPa = 133.32

    # From El Bouri + Payne, we know
    #Kc = 4.28e-4 # in mm^3 s kg^-1
    Ka = 1.234
    Kv = 2.468

    Kc = 3.3e-3 # from smith et al (2014) "transmural variation and anisotropy..."

    # permeability and viscosity
    mu_blood = 2.67e-3 #Pa.s
    kc = Kc*1e-9*mu_blood
    ka = Ka*1e-9*mu_blood
    kv = Kv*1e-9*mu_blood

    # Resitance coeffs from Vinje (2020)
    mlTomm = 1000.0
    factor_resistance = mmHgToPa*60.0/mlTomm

    # Resistances in perivascular spaces
    R_pa = 1.14 * factor_resistance
    R_pv = 1.75e-3 * factor_resistance
    R_pc = 32.24 * factor_resistance

    # Resistance across membranes
    R_pa_e = 0.57 * factor_resistance
    R_pv_e = 0.64 * factor_resistance
    R_c_pc = 125.31 * factor_resistance

    ### Permeabilities
    # 1st step determine the ration
    # We know the permeability in the ECS and its resistance
    kappa_ECS = 2.0e-11 # mm^2
    R_ECS = 0.57 * factor_resistance
    mu_CSF = 7.e-4 # Pa.s
    ratio_L_A = kappa_ECS*R_ECS/mu_CSF
    print('ratio L/A = ' +str(ratio_L_A))

    kappa_a = ka*1.e6
    kappa_v = kv*1.e6
    kappa_c = kc*1.e6

    kappa_pa = mu_CSF/R_pa*ratio_L_A
    kappa_pv = mu_CSF/R_pv*ratio_L_A
    kappa_pc = mu_CSF/R_pc*ratio_L_A

    print("###  Permeabilities  ###\n")
    print("kappa_a = " +str(kappa_a))
    print("kappa_v = " +str(kappa_v))
    print("kappa_c = " +str(kappa_c))

    print("kappa_pa = " +str(kappa_pa))
    print("kappa_pv = " +str(kappa_pv))
    print("kappa_pc = " +str(kappa_pc))

    print("kappa_ECS = " +str(kappa_ECS))

    ## Transfer coefficients
    # across BBB
    # Hydraulic conductivity
    L_a_e = 9.1e-10
    L_c_e = 1.0e-10
    L_v_e = 2.0e-11

    # surface to volume ration of vessels
    surface_v_to_volume = 3.  #mm^-1
    surface_a_to_volume = 3. #mm^-1
    surface_c_to_volume = 9.  #mm^-1

    # compute transfer from blood vessels to ECS
    gamma_a_e = L_a_e*surface_a_to_volume
    gamma_v_e = L_v_e*surface_v_to_volume
    gamma_c_e = L_c_e*surface_c_to_volume



    print("\n###  Transfer  ###\n")
    print("gamma_a_e = " +str(gamma_a_e))
    print("gamma_v_e = " +str(gamma_v_e))
    print("gamma_c_e = " +str(gamma_c_e))

    # PVS to tissue
    gamma_pa_e = 1.0/(R_pa_e*volume_brain_human)
    gamma_pv_e = 1.0/(R_pv_e*volume_brain_human)
    R_pc_e = (1/(gamma_c_e*volume_brain_human) - R_c_pc) # compute resistance
    gamma_pc_e = 1.0/(R_pc_e*volume_brain_human)

    print("gamma_pa_e = " +str(gamma_pa_e))
    print("gamma_pv_e = " +str(gamma_pv_e))
    print("gamma_pc_e = " +str(gamma_pc_e))

    # From blood to PVS
    # Compute resistance from blood vessel to corresponding PVS
    R_a_pa = 1/(gamma_a_e*volume_brain_human)-R_pa_e
    R_v_pv = 1/(gamma_v_e*volume_brain_human)-R_pv_e

    gamma_a_pa = 1.0/(R_a_pa*volume_brain_human)
    gamma_v_pv = 1.0/(R_v_pv*volume_brain_human)
    gamma_c_pc = 1.0/(R_c_pc*volume_brain_human)

    print("gamma_a_pa = " +str(gamma_a_pa))
    print("gamma_v_pv = " +str(gamma_v_pv))
    print("gamma_c_pc = " +str(gamma_c_pc))

    ## Transfer along vessels
    # From PVS arteries to PVS capillaries and PVS capillaries to PVS veins
    deltap_pa_pc = 1.0*mmHgToPa # difference of pressure between PVS arteries and PVS capillaries
    deltap_pc_pv = 0.25*mmHgToPa # difference of pressure between PVS arteries and PVS capillaries
    B_PVS = 3.38/60 # (mm^3/s) flow in PVS

    # transfer coefficients
    gamma_pa_pc = B_PVS/(deltap_pa_pc*volume_brain_rat)
    gamma_pc_pv = B_PVS/(deltap_pc_pv*volume_brain_rat)


    # From arteries to capillaries and capillaries to veins
    deltap_a_c = 40.0*mmHgToPa # difference of pressure between arteries and capillaries
    deltap_c_v = 13.0*mmHgToPa # difference of pressure between capillaries and veins
    B_blood = 116.*1.e3/100*2./60 # flow in blood network

    gamma_a_c = B_blood/(deltap_a_c*volume_brain_rat)
    gamma_c_v = B_blood/(deltap_c_v*volume_brain_rat)



    print("\ngamma_a_c = " + str(gamma_a_c))
    print("gamma_c_v = " +str(gamma_c_v))
    print("gamma_pa_pc = " +str(gamma_pa_pc))
    print("gamma_pc_pv = " +str(gamma_pc_pv))

    # Coefficient boundaries
    print("\n### Hydraulic conductivities at pial surface")
    surface_brain = 1750*1.e2
    gamma_e_SAS = 1/(2*R_pa)

    print("L_e_SAS = " + str(gamma_e_SAS/surface_brain))
    print("L_pa_SAS = " + str(gamma_pa_e/surface_brain*volume_brain_human))
    L_e_SAS = gamma_e_SAS/surface_brain
    L_pa_SAS = gamma_pa_e/surface_brain*volume_brain_human


    ### Diffusive Permeabilities
    print("\n### Computing diffusive permeability coefficients ###")
    print("## For Inulin ##")
    print("# Across the Astrocyte endfeet barrier# \n")
    # Computing diffusive permeabilities from Fu article  "A model for the blood–brain barrier permeability to water and small solutes"

    ## For Inulin
    # Computing effective diffusion in AEF pores at the capillary level
    Dv = 10000.0e-6 # diameter capillary
    D_Free_Inulin = 2.98e-4 # Free diffusion coefficient # FIXED: paper reports 2.98e-6 cm^2/s -> 2.98e-4 mm^2/s
    a_Inulin = 15.2e-7 # mm

    L_AEF = 1000.0e-6 # Width of AEF (assumed constant at the different levels)
    B_AEF = 10.0e-6 # mm can inccrease up to 1000e-6 (We assume very small because we are at the capillary level)
    beta = a_Inulin/B_AEF
    D_eff_AEF_Inulin = D_Free_Inulin*(1-2.10444*beta +2.08877*pow(beta,3) - 0.094813*pow(beta,5) - 1.372*pow(beta,6)) # FIXED: removed errouneous beta**6 in second term.

    R_AEF_Inulin = L_AEF/(2*B_AEF*D_eff_AEF_Inulin)

    P_AEF_Inulin = 1/R_AEF_Inulin*1/(np.pi*Dv)
    lambda_pc_e_Inulin = P_AEF_Inulin*surface_c_to_volume

    print("P_pc,e_Inulin = "  +str(P_AEF_Inulin))


    # Computing effective diffusion in AEF pores at the venous and arterial levels
    Dv = 50000.0e-6 # diameter arteriole or venule
    Dv = 38.e-3

    B_AEF = 250.0e-6 # mm can decrease to 2.5e-6 (We assume large because we are at the venous and arterial levels)
    beta = a_Inulin/B_AEF
    D_eff_AEF_Inulin = D_Free_Inulin*(1-2.10444*beta +2.08877*pow(beta,3) - 0.094813*pow(beta,5) - 1.372*pow(beta,6)) # FIXED: removed errouneous beta**6 in second term.

    R_AEF_Inulin = L_AEF/(2*B_AEF*D_eff_AEF_Inulin)

    P_AEF_Inulin = 1/R_AEF_Inulin*1/(np.pi*Dv)
    lambda_pa_e_Inulin = P_AEF_Inulin*surface_a_to_volume
    lambda_pv_e_Inulin = P_AEF_Inulin*surface_v_to_volume
    print("P_pa,e_Inulin = P_pv,e_Inulin ="  +str(P_AEF_Inulin))

    print("\nlambda_pa,e_Inulin = "  +str(lambda_pa_e_Inulin))
    print("lambda_pc,e_Inulin = "  +str(lambda_pc_e_Inulin))
    print("lambda_pv,e_Inulin = "  +str(lambda_pv_e_Inulin))

    print("\n## For Amyloid-beta ##")
    print("# Across the BBB#") # We need to bring down all the different layers
    ## Capillary level
    # Astrocyte process
    L_AEF = 1000.0e-6 # Width of AEF (assumed constant at the different levels)
    B_AEF = 10.0e-6 # mm can inccrease up to 1000e-6 (We assume very small because we are at the capillary level)
    a_AB = 0.9e-6 # mm
    beta = a_AB/B_AEF
    D_Free_AB = 1.8e-6
    D_eff_AEF_AB = D_Free_AB*(1-2.10444*pow(beta,6) +2.08877*pow(beta,3) - 0.094813*pow(beta,5) - 1.372*pow(beta,6))


    R_AEF_AB = L_AEF/(2*B_AEF*D_eff_AEF_AB)
    P_AEF_AB = 1/R_AEF_AB*1/(np.pi*Dv)

    # We already know
    lambda_pc_e_AB = P_AEF_AB*surface_c_to_volume

    # SGL layer
    B= 9.0e-6 # width of inter-endothelial cleft
    L_f = 250.0e-6 # Glycocalyx layer thickness
    r_f = 6e-6 # radius of fibers
    eps = 1-0.326 # volume fraction of void in glycocalyx layer
    D_f_AB =  D_Free_AB*(np.exp( -pow(1-eps,0.5)*(1+ a_AB/(r_f )) )) # unordered fibers see Michel 1999

    R_f = L_f/(2*B*D_f_AB)

    # inter-endothelial cleft
    L_1 = 350.0e-6 # Distance between the tight junction strand and the front of the inter-endothelial cleft

    beta = a_AB/B
    D_c_AB =  D_Free_AB*(1-2.10444*pow(beta,6) +2.08877*pow(beta,3) - 0.094813*pow(beta,5) - 1.372*pow(beta,6))
    R_1_AB = L_1/(2*B*D_c_AB)

    # small silt in the junction strand
    L_jun = 11.0e-6 # Thickness of tight junction strand
    B_s = 2.5e-6  # Width of the small slit of the tight junction strand (Allt and Lawrenson, 1997; Cassella et al., 1997)
    beta = a_AB/B_s
    D_sl_AB =  D_Free_AB*(1-2.10444*pow(beta,6) +2.08877*pow(beta,3) - 0.094813*pow(beta,5) - 1.372*pow(beta,6))
    R_2_AB = L_jun/(2*B_s*D_sl_AB)

    # small slit in the junction strand
    L = 700.0e-6 # Total length of the cleft region (Schulze and Firth, 1992)
    R_3_AB = (L-L_1-L_f)/(2*B*D_c_AB)

    # Basement membrane (fibrous membrane)
    W_a = 2500.0e-6 # Length of astrocyte foot processes (Farkas and Luiten, 2001)
    r_b = 2.7e-6 # radius of BM fibers Michel (1999)
    eps_b = 1-0.5 # void volume fraction in BM Michel (1999)
    L_b = 30.0e-6 # Thickness of the basement membrane (Paulson and Newman, 1987; Farkas and Luiten, 2001)
    D_b_AB =  D_Free_AB*(np.exp( -pow(1-eps_b,0.5)*(1+ a_AB/(r_b )) )) # unordered fibers see Michel 1999
    R_4_AB = W_a/(4*L_b*D_b_AB)

    # Permeability across total BBB at the capillary level
    Dv = 10.0e-3
    P_BBB_AB = 1/(R_AEF_AB + R_f  + R_1_AB + R_2_AB + R_3_AB + R_4_AB) * 1/(np.pi * Dv)


    lambda_c_e_AB = P_BBB_AB * surface_c_to_volume
    P_BBB_del_AEF_AB  = 1/(R_f  + R_1_AB + R_2_AB + R_3_AB + R_4_AB) * 1/(np.pi * Dv)
    lambda_c_pc_AB = P_BBB_del_AEF_AB * surface_c_to_volume


    ### Now the arterial space, we assume that the junctions are tight because the BBB is known to be less permeable at this level
    # Astrocyte process
    L_AEF = 1000.0e-6 # Width of AEF (assumed constant at the different levels)
    B_AEF = 500.0e-6 # mm can inccrease up to 1000e-6 (We assume large because we are at the artery level)
    beta = a_AB/B_AEF
    D_eff_AEF_AB = D_Free_AB*(1-2.10444*pow(beta,6) +2.08877*pow(beta,3) - 0.094813*pow(beta,5) - 1.372*pow(beta,6))

    R_AEF_AB = L_AEF/(2*B_AEF*D_eff_AEF_AB)
    P_AEF_AB = 1/R_AEF_AB*1/(np.pi*Dv)

    # We already know
    lambda_pa_e_AB = P_AEF_AB*surface_a_to_volume

    # SGL layer
    B= 9.0e-6 # width of inter-endothelial cleft
    L_f = 400.0e-6 # Glycocalyx layer thickness (assumed to be large at this level)
    r_f = 6e-6 # radius of fibers
    eps = 1-0.326 # volume fraction of void in glycocalyx layer
    D_f_AB =  D_Free_AB*(np.exp( -pow(1-eps,0.5)*(1+ a_AB/(r_f )) )) # unordered fibers see Michel 1999

    R_f = L_f/(2*B*D_f_AB)

    # inter-endothelial cleft
    L_1 = 350.0e-6 # Distance between the tight junction strand and the front of the inter-endothelial cleft

    beta = a_AB/B
    D_c_AB =  D_Free_AB*(1-2.10444*pow(beta,6) +2.08877*pow(beta,3) - 0.094813*pow(beta,5) - 1.372*pow(beta,6))
    R_1_AB = L_1/(2*B*D_c_AB)

    # small silt in the junction strand
    L_jun = 11.0e-6 # Thickness of tight junction strand
    B_s = 0.5e-6  # Width of the small slit of the tight junction strand (Allt and Lawrenson, 1997; Cassella et al., 1997)
    beta = a_AB/B_s
    D_sl_AB =  D_Free_AB*(1-2.10444*pow(beta,6) +2.08877*pow(beta,3) - 0.094813*pow(beta,5) - 1.372*pow(beta,6))
    R_2_AB = L_jun/(2*B_s*D_sl_AB)

    # small slit in the junction strand
    L = 700.0e-6 # Total length of the cleft region (Schulze and Firth, 1992)
    R_3_AB = (L-L_1-L_f)/(2*B*D_c_AB)

    # Basement membrane (fibrous membrane)
    W_a = 2500.0e-6 # Length of astrocyte foot processes (Farkas and Luiten, 2001)
    r_b = 2.7e-6 # radius of BM fibers Michel (1999)
    eps_b = 1-0.5 # void volume fraction in BM Michel (1999)
    L_b = 80.0e-6 # Thickness of the basement membrane (Paulson and Newman, 1987; Farkas and Luiten, 2001)
    D_b_AB =  D_Free_AB*(np.exp( -pow(1-eps_b,0.5)*(1+ a_AB/(r_b )) )) # unordered fibers see Michel 1999
    R_4_AB = W_a/(4*L_b*D_b_AB)

    # Permeability across total BBB at the capillary level
    Dv = 50.0e-3
    P_BBB_AB = 1/(R_AEF_AB + R_f  + R_1_AB + R_2_AB + R_3_AB + R_4_AB) * 1/(np.pi * Dv)

    lambda_a_e_AB = P_BBB_AB * surface_a_to_volume

    P_BBB_del_AEF_AB = 1/( R_f  + R_1_AB + R_2_AB + R_3_AB + R_4_AB) * 1/(np.pi * Dv)

    lambda_a_pa_AB = P_BBB_del_AEF_AB * surface_a_to_volume

    ### Now the VENOUS space, we assume that the junctions are loose because the BBB is known to be more permeable at this level (less claudin expression)
    # Astrocyte process
    L_AEF = 1000.0e-6 # Width of AEF (assumed constant at the different levels)
    B_AEF = 500.0e-6 # mm can inccrease up to 1000e-6 (We assume large because we are at the artery level)
    beta = a_AB/B_AEF
    D_eff_AEF_AB = D_Free_AB*(1-2.10444*pow(beta,6) +2.08877*pow(beta,3) - 0.094813*pow(beta,5) - 1.372*pow(beta,6))

    R_AEF_AB = L_AEF/(2*B_AEF*D_eff_AEF_AB)
    P_AEF_AB = 1/R_AEF_AB*1/(np.pi*Dv)

    # We already know
    lambda_pv_e_AB = P_AEF_AB*surface_v_to_volume

    # SGL layer
    B= 9.0e-6 # width of inter-endothelial cleft
    L_f = 100.0e-6 # Glycocalyx layer thickness (assumed to be large at this level)
    r_f = 6e-6 # radius of fibers
    eps = 1-0.326 # volume fraction of void in glycocalyx layer
    D_f_AB =  D_Free_AB*(np.exp( -pow(1-eps,0.5)*(1+ a_AB/(r_f )) )) # unordered fibers see Michel 1999

    R_f = L_f/(2*B*D_f_AB)

    # inter-endothelial cleft
    L_1 = 350.0e-6 # Distance between the tight junction strand and the front of the inter-endothelial cleft

    beta = a_AB/B
    D_c_AB =  D_Free_AB*(1-2.10444*pow(beta,6) +2.08877*pow(beta,3) - 0.094813*pow(beta,5) - 1.372*pow(beta,6))
    R_1_AB = L_1/(2*B*D_c_AB)

    # small silt in the junction strand
    L_jun = 11.0e-6 # Thickness of tight junction strand
    B_s = 10.0e-6  # Width of the small slit of the tight junction strand (Allt and Lawrenson, 1997; Cassella et al., 1997)
    beta = a_AB/B_s
    D_sl_AB =  D_Free_AB*(1-2.10444*pow(beta,6) +2.08877*pow(beta,3) - 0.094813*pow(beta,5) - 1.372*pow(beta,6))
    R_2_AB = L_jun/(2*B_s*D_sl_AB)

    # small slit in the junction strand
    L = 700.0e-6 # Total length of the cleft region (Schulze and Firth, 1992)
    R_3_AB = (L-L_1-L_f)/(2*B*D_c_AB)

    # Basement membrane (fibrous membrane)
    W_a = 2500.0e-6 # Length of astrocyte foot processes (Farkas and Luiten, 2001)
    r_b = 2.7e-6 # radius of BM fibers Michel (1999)
    eps_b = 1-0.5 # void volume fraction in BM Michel (1999)
    L_b = 20.0e-6 # Thickness of the basement membrane (Paulson and Newman, 1987; Farkas and Luiten, 2001)
    D_b_AB =  D_Free_AB*(np.exp( -pow(1-eps_b,0.5)*(1+ a_AB/(r_b )) )) # unordered fibers see Michel 1999
    R_4_AB = W_a/(4*L_b*D_b_AB)

    # Permeability across total BBB at the capillary level
    Dv = 50.0e-3
    P_BBB_AB = 1/(R_AEF_AB + R_f  + R_1_AB + R_2_AB + R_3_AB + R_4_AB) * 1/(np.pi * Dv)

    lambda_v_e_AB = P_BBB_AB * surface_v_to_volume
    P_BBB_del_AEF_AB = 1/(R_f  + R_1_AB + R_2_AB + R_3_AB + R_4_AB) * 1/(np.pi * Dv)
    lambda_v_pv_AB = P_BBB_del_AEF_AB * surface_v_to_volume




    print("\nlambda_a_e_AB = "+str(lambda_a_e_AB))
    print("lambda_v_e_AB = "+str(lambda_v_e_AB))
    print("lambda_c_e_AB = "+str(lambda_c_e_AB))

    print("\n## Across the AEF")

    print("\nlambda_pa_e_AB = "+str(lambda_pa_e_AB))
    print("lambda_pv_e_AB = "+str(lambda_pv_e_AB))
    print("lambda_pc_e_AB = "+str(lambda_pc_e_AB))

    print("\n## Across the BBB without the AEF")
    print("\nlambda_a_pa_AB = "+str(lambda_a_pa_AB))
    print("lambda_v_pv_AB = "+str(lambda_v_pv_AB))
    print("lambda_c_pc_AB = "+str(lambda_c_pc_AB))

    ### Convective mass exchange
    ## Note The reflection coefficient is often found from microvessels but it depends on the membrane itself.
    ## Whether this value changes from the AEF membrane only is unclear. Might need further investigation.
    ## Inulin
    sigma_reflect = 0.2


    tilde_gamma_pa_e_Inulin = gamma_pa_e*(1-sigma_reflect)
    tilde_gamma_pv_e_Inulin = gamma_pv_e*(1-sigma_reflect)
    tilde_gamma_pc_e_Inulin = gamma_pc_e*(1-sigma_reflect)


    # For connected channels, we assume that there is not a membrane (so no diffusive transfer)
    sigma_reflect = 0.0

    tilde_gamma_pa_pc_Inulin = gamma_pa_pc
    tilde_gamma_pc_pv_Inulin = gamma_pc_pv

    ## Amyloid beta
    sigma_reflect = 0.1

    tilde_gamma_a_e_AB = gamma_a_e*(1-sigma_reflect)
    tilde_gamma_v_e_AB = gamma_v_e*(1-sigma_reflect)
    tilde_gamma_c_e_AB = gamma_c_e*(1-sigma_reflect)

    tilde_gamma_pa_e_AB = gamma_pa_e*(1-sigma_reflect)
    tilde_gamma_pv_e_AB = gamma_pv_e*(1-sigma_reflect)
    tilde_gamma_pc_e_AB = gamma_pc_e*(1-sigma_reflect)

    tilde_gamma_a_pa_AB = gamma_a_pa*(1-sigma_reflect)
    tilde_gamma_v_pv_AB = gamma_v_pv*(1-sigma_reflect)
    tilde_gamma_c_pc_AB = gamma_c_pc*(1-sigma_reflect)


    # For connected channels, we assume that there is not a membrane (so no diffusive transfer)
    sigma_reflect = 0.0

    tilde_gamma_a_c_AB = gamma_a_c
    tilde_gamma_c_v_AB = gamma_c_v
    tilde_gamma_pa_pc_AB = gamma_pa_pc
    tilde_gamma_pc_pv_AB = gamma_pc_pv

    print("\n### Convective mass transfer ###\n")
    print("## For Inulin\n")
    print("tilde_gamma_pa_e_Inulin = "+str(tilde_gamma_pa_e_Inulin))
    print("tilde_gamma_pv_e_Inulin = "+str(tilde_gamma_pv_e_Inulin))
    print("tilde_gamma_pc_e_Inulin = "+str(tilde_gamma_pc_e_Inulin))

    print("\ntilde_gamma_pa_pc_Inulin = "+str(tilde_gamma_pa_pc_Inulin))
    print("tilde_gamma_pc_pv_Inulin = "+str(tilde_gamma_pc_pv_Inulin))




    print("\n## For Amyloid-Beta\n")

    print("tilde_gamma_a_e_AB = " + str(tilde_gamma_a_e_AB))
    print("tilde_gamma_v_e_AB = " + str(tilde_gamma_v_e_AB))
    print("tilde_gamma_c_e_AB  = "+ str(tilde_gamma_c_e_AB))

    print("tilde_gamma_pa_e_AB = " + str(tilde_gamma_pa_e_AB))
    print("tilde_gamma_pv_e_AB = " + str(tilde_gamma_pv_e_AB))
    print("tilde_gamma_pc_e_AB = " + str(tilde_gamma_pc_e_AB))

    print("tilde_gamma_a_pa_AB = " + str(tilde_gamma_a_pa_AB))
    print("tilde_gamma_v_pv_AB = " + str(tilde_gamma_v_pv_AB))
    print("tilde_gamma_c_pc_AB = " + str(tilde_gamma_c_pc_AB))


    # Save to dictionary

    volume_fraction_of_blood = 3.29/100.
    arteries_frac = 0.21
    veins_frac = 0.46
    capillaries_frac = 0.33
    PVS_fraction = 1./100.

    main_dict["phi_e"] = 0.14
    main_dict["phi_a"] = volume_fraction_of_blood*arteries_frac
    main_dict["phi_v"] = volume_fraction_of_blood*veins_frac
    main_dict["phi_c"] = volume_fraction_of_blood*capillaries_frac
    main_dict["phi_pa"] = PVS_fraction*arteries_frac
    main_dict["phi_pv"] = PVS_fraction*veins_frac
    main_dict["phi_pc"] = PVS_fraction*capillaries_frac

    print("\n### Baseline porosities ###\n")
    print("phi_e = "+str(main_dict["phi_e"]))
    print("phi_a = "+str(main_dict["phi_a"]))
    print("phi_v = "+str(main_dict["phi_v"]))
    print("phi_c = "+str(main_dict["phi_c"]))
    print("phi_pa = "+str(main_dict["phi_pa"]))
    print("phi_pv = "+str(main_dict["phi_pv"]))
    print("phi_pc = "+str(main_dict["phi_pc"]))


    main_dict["mu_blood"] = mu_blood
    main_dict["mu_csf"] = mu_CSF

    main_dict["gamma_a_e"] = gamma_a_e
    main_dict["gamma_v_e"] = gamma_v_e
    main_dict["gamma_c_e"] = gamma_c_e
    main_dict["gamma_pa_e"] = gamma_pa_e
    main_dict["gamma_pv_e"] = gamma_pv_e
    main_dict["gamma_pc_e"] = gamma_pc_e
    main_dict["gamma_a_pa"] = gamma_a_pa
    main_dict["gamma_v_pv"] = gamma_v_pv
    main_dict["gamma_c_pc"] = gamma_c_pc
    main_dict["gamma_a_c"] = gamma_a_c
    main_dict["gamma_c_v"] = gamma_c_v
    main_dict["gamma_pa_pc"] = gamma_pa_pc
    main_dict["gamma_pc_pv"] = gamma_pc_pv
    main_dict["kappa_a"] = kappa_a
    main_dict["kappa_v"] = kappa_v
    main_dict["kappa_c"] = kappa_c
    main_dict["kappa_pa"] = kappa_pa
    main_dict["kappa_pv"] = kappa_pv
    main_dict["kappa_pc"] = kappa_pc
    main_dict["kappa_e"] = kappa_ECS
    main_dict["lambda_pa_e_Inulin"] = lambda_pa_e_Inulin
    main_dict["lambda_pv_e_Inulin"] = lambda_pv_e_Inulin
    main_dict["lambda_pc_e_Inulin"] = lambda_pc_e_Inulin
    main_dict["tilde_gamma_pa_e_Inulin"] = tilde_gamma_pa_e_Inulin
    main_dict["tilde_gamma_pc_e_Inulin"] = tilde_gamma_pc_e_Inulin
    main_dict["tilde_gamma_pv_e_Inulin"] = tilde_gamma_pv_e_Inulin
    main_dict["tilde_gamma_pa_pc_Inulin"] = tilde_gamma_pa_pc_Inulin
    main_dict["tilde_gamma_pc_pv_Inulin"] = tilde_gamma_pc_pv_Inulin

    main_dict["L_e_SAS"] = L_e_SAS
    main_dict["L_pa_SAS"] = L_pa_SAS

    main_dict["blood_flow"] = B_blood


    if type_exp == "ECS enhanced":
        main_dict["phi_e"] = 0.23
        main_dict["kappa_e"] = kappa_ECS*5.5
    elif type_exp == "ECS+PVS enhanced":
        main_dict["phi_e"] = 0.23
        main_dict["kappa_e"] = kappa_ECS*5.5
        main_dict["phi_pa"] = main_dict["phi_pa"]*4.0
        main_dict["phi_pv"] = main_dict["phi_pv"]*4.0
        main_dict["phi_pc"] = main_dict["phi_pc"]*4.0
        main_dict["kappa_pa"] = kappa_pa*16.
        main_dict["kappa_pv"] = kappa_pv*16.
        main_dict["kappa_pc"] = kappa_pc *16.
    elif type_exp == "PVS only enhanced":
        main_dict["phi_pa"] = main_dict["phi_pa"]*4.0
        main_dict["phi_pv"] = main_dict["phi_pv"]*4.0
        main_dict["phi_pc"] = main_dict["phi_pc"]*4.0
        main_dict["kappa_pa"] = kappa_pa*16.
        main_dict["kappa_pv"] = kappa_pv*16.
        main_dict["kappa_pc"] = kappa_pc *16.

    return main_dict

if __name__ == "__main__":
    compute_coefficients(None)
