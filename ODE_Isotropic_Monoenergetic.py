from Auxiliar_Functions import *
from scipy.integrate import odeint

def func_Isotropic_Monoenergetic(y, time, params):
    omega,mu_opt,mu_0,lamb_opt,lamb_0,n_dim= params  # unpack parameters
    
    B=np.array(B_vec(n_dim,theta_31))
    L=np.array(L_vec(n_dim))
    
    r=time/from_eV_to_1_over_km #From eV⁻¹ to km
    mu=mu_supernova(r,mu_opt,mu_0)
    lamb=lambda_supernova(r,lamb_opt,lamb_0)

    derivs=[]
    P_aux= cross_prod(y[0:n_dim],(B*omega+L*lamb+mu*y[n_dim:]))
    for k in range(n_dim):
        derivs.append(P_aux[k])
    
    #nu_bar 
    P_aux= cross_prod(y[n_dim:],(-1*B*omega+L*lamb-mu*y[0:n_dim]))
    for k in range(n_dim):
        derivs.append(P_aux[k])
        
    return derivs


def solver_Isotropic_Monoenergetic(P,E,r_i,r_f,mass_ord,mu_opt,mu_0,lamb_opt="no",lamb_0=0,n_f=2):

    omega=delta_m2_31/(2*E*10**6) #eV
    r_step = (2*np.pi/max(omega,mu_0))/200 #eV⁻¹
    r_i = r_i*from_eV_to_1_over_km
    r_f = r_f*from_eV_to_1_over_km
    r = np.arange(r_i,r_f,r_step) #eV⁻¹
    n_dim=(n_f**2)-1
    
    if mass_ord=="NH": 
        params=omega,mu_opt,mu_0,lamb_opt,lamb_0,n_dim
    elif mass_ord=="IH":
        omega=-1*omega
        params=omega,mu_opt,mu_0,lamb_opt,lamb_0,n_dim
    else:
        print("Not a mass ordering option!")
        return 0

    psoln  = odeint(func_Isotropic_Monoenergetic, P, r, args=(params,))
    psoln_trans=np.transpose(psoln)
    P_nu=psoln_trans[0:3]
    P_nubar=psoln_trans[3:6]
    
    H_vac,H_nue=[],[]
    r=r/from_eV_to_1_over_km #From eV⁻¹ to km
    
    for r_i in r:
        B=np.array(B_vec(n_dim,theta_31))
        L=np.array(L_vec(n_dim))
        lamb=lambda_supernova(r,lamb_opt,lamb_0)
        H_vac.append(omega*B)
        H_nue.append(lamb*L)
        
    H_vac=np.transpose(H_vac)  
    H_nue=np.transpose(H_nue)   
    
    return P_nu,P_nubar,H_vac,H_nue,r