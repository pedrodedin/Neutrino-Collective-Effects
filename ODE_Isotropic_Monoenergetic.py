from Auxiliar_Functions import *
from scipy.integrate import odeint

def func_Isotropic_Monoenergetic(y, time, params):
#     nu1, nu2, nu3, antinu1, antinu2, antinu3 = y      # unpack current values of y
    omega,mu_opt,mu_0,lamb_opt,lamb_0,n_dim= params  # unpack parameters
    
    #B1,B2,B3=omega*np.array(B_vec(n_dim,theta_31))
    B=np.array(B_vec(n_dim,theta_31))
    L=np.array(L_vec(n_dim))
    
    r=time/from_eV_to_1_over_km #From eV⁻¹ to km
    mu=mu_supernova(r,mu_opt,mu_0)
    lamb=lambda_supernova(r,lamb_opt,lamb_0)

#     derivs = [nu2*B3-nu3*B2 +mu*(nu2*antinu3-nu3*antinu2) ,      # list of dy/dt=f functions
#               nu3*B1-nu1*B3 +mu*(nu3*antinu1-nu1*antinu3) ,
#               nu1*B2-nu2*B1 +mu*(nu1*antinu2-nu2*antinu1) ,
#               -1*(antinu2*B3-antinu3*B2) +mu*(nu2*antinu3-nu3*antinu2),      
#               -1*(antinu3*B1-antinu1*B3) +mu*(nu3*antinu1-nu1*antinu3),
#               -1*(antinu1*B2-antinu2*B1) +mu*(nu1*antinu2-nu2*antinu1)] 
    #nu 
    #print(y[0:n_dim],y[n_dim:])
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
        params=-1*omega,mu_opt,mu_0,lamb_opt,lamb_0,n_dim
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