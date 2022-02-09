from Auxiliar_Functions import *
from scipy.integrate import odeint
R_sol = 6.96340*10**8 #m
R_sol= R_sol*(8*10**5) #eV⁻¹
    
def N_e_sol(r):
    r_0=R_sol/10.54
    N_0=245 #N_A cm⁻³
    N_0=N_A*N_0*(1.973*10**(-5))**3#eV³ 
    N = N_0*math.exp(-r/r_0)
    
    return N

#EDO
def f_solar(y, t, params):
    S1, S2, S3 = y      # unpack current values of y
    E= params  # unpack parameters
    
    N_e=N_e_sol(t)
    #Mixing parameters in matter
    theta_M=theta_eff(delta_m2_21,theta_21,Acc(N_e,E))
    delta_m2_M=delta_m2_eff(delta_m2_21,theta_21,Acc(N_e,E))
    omega_M=delta_m2_M/(2*E)
    
    B1,B2,B3 = omega_M*np.array(B_vec(3,theta_M))
    derivs = [S2*B3-S3*B2,      # list of dy/dt=f functions
              S3*B1-S1*B3,
              S1*B2-S2*B1]
    return derivs


def solver_solar(E,r_i,r_f):
    #E [MeV], r_i [R_solar],  r_i [R_solar]
    E=E*10**6 #eV
    #omega_0=delta_m2_21/(2*E*10**6)
    r_step=0.000001 # units of R_solar
    r_per_R_sol=np.arange(r_i,r_f,r_step) # units of R_solar
    r=R_sol*r_per_R_sol #eV⁻¹

    # Initial values
    S1_0 = 0   
    S2_0 = 0     
    S3_0 = 1.0  #Electronic Neutrino

    # Initial conditions for ODE solver
    y0 = [S1_0, S2_0, S3_0]

    # Call the ODE solver
    psoln=odeint(f_solar, y0, r, args=(E,))    
    psoln_trans=np.transpose(psoln,(1,0))
    
    S1=psoln_trans[0]
    S2=psoln_trans[1]
    S3=psoln_trans[2]
    S=np.array([S1,S2,S3])
    
    B=[]
    for r_i in r:
        N_e=N_e_sol(r_i)
        #Mixing parameters in matter
        theta_M=theta_eff(delta_m2_21,theta_21,Acc(N_e,E))
        delta_m2_M=delta_m2_eff(delta_m2_21,theta_21,Acc(N_e,E))
        omega_M=delta_m2_M/(2*E)
        B.append(omega_M*np.array(B_vec(3,theta_M)))
    B=np.transpose(B,(1,0))
    
    return S,B,r_per_R_sol