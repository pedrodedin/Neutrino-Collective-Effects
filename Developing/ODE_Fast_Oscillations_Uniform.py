from Auxiliar_Functions import *
from scipy.integrate import odeint,solve_ivp

##### Angular Distributions ######
def angular_dist_Ian(theta,a,b):
    #Neutrino Flavor Pendulum Reloaded: The Case of Fast Pairwise Conversion - Padilla-Gay, Irene Tamborra, Georg G. Raffelt 
    # https://arxiv.org/abs/2109.14627
    return 0.45 - a +(0.1/b)*np.exp(-1*(1-np.cos(theta))**2/(2*b**2))

def theta_dist(theta,nu_type,case):
    if nu_type=="nu_x" or nu_type=="nu_x_bar":
        return np.where(theta<-100,0.0,0.0)
    
    # New Developments in Flavor Evolution of a Dense Neutrino Gas - Irene Tamborra, Shashank Shalgar
    #https://arxiv.org/abs/2011.01948
    if case==1:
        if nu_type=="nu_e":
            return np.where(theta<-100,0,0.5)
        if nu_type=="nu_e_bar":
            return np.where(theta<math.pi/3,1,0.25)
    elif case==2:
        if nu_type=="nu_e":
            return np.where(theta<-100,0,0.5)
        if nu_type=="nu_e_bar":
            return np.where(theta<math.pi/4,1,0.5)
    elif case==3:
        if nu_type=="nu_e":
            return np.where(theta<-100,0,0.5)
        if nu_type=="nu_e_bar":
            return np.where(theta<math.pi/3,0.5*1.5,0.5)
        
    #Neutrino Flavor Pendulum Reloaded: The Case of Fast Pairwise Conversion - Padilla-Gay, Irene Tamborra, Georg G. Raffelt 
    # https://arxiv.org/abs/2109.14627
    elif case=='A':
        if nu_type=="nu_e":
            return np.where(theta<-100,0,0.5)
        if nu_type=="nu_e_bar":
            return angular_dist_Ian(theta,0,0.4)
    elif case=='B':
        if nu_type=="nu_e":
            return np.where(theta<-100,0,0.5)
        if nu_type=="nu_e_bar":
            return angular_dist_Ian(theta,0.02,0.4)       
    elif case=='C':
        if nu_type=="nu_e":
            return np.where(theta<-100,0,0.5)
        if nu_type=="nu_e_bar":
            return angular_dist_Ian(theta,0.02,0.6)
    elif case=='D':
        if nu_type=="nu_e":
            return np.where(theta<-100,0,0.5)
        if nu_type=="nu_e_bar":
            return angular_dist_Ian(theta,0.06,0.2)
    
    else:
        print("Not a valid angular distribution!")

        
####################### ODE ################################
####### Initiate #######
def initiate_FFC_uniform(nu_types,t_i,t_f,E_nu,delta_m2,mu_0,theta_bins,case,seed=0):
    y0=[] #Initial state
    flavor_sign=[1,-1]
    omega=delta_m2/(2*E_nu*10**6) #eV 
    rho_ex=0
    if delta_m2==0:
        rho_ex=seed

    cos_theta_i,cos_theta_f=-1,1
    cos_theta_vec=np.linspace(cos_theta_i,cos_theta_f,theta_bins)
    cos_theta_step=(cos_theta_f-cos_theta_i)/theta_bins
    n_theta=len(cos_theta_vec)
    n_f=len(nu_types)
    n_dim=(n_f**2)-1
    n_antipart=2

    index_order=[n_theta,n_f,n_dim,n_antipart]
    
    for i in range(n_theta):        
      for j in range(n_f):
        #nu
        nu_spec=theta_dist(np.arccos(cos_theta_vec[i]),nu_types[j],case)*cos_theta_step
        y0.append(rho_ex)
        y0.append(rho_ex)
        y0.append(flavor_sign[j]*nu_spec)
        #nubar
        nu_spec=theta_dist(np.arccos(cos_theta_vec[i]),nu_types[j]+'_bar',case)*cos_theta_step
        y0.append(rho_ex)
        y0.append(rho_ex)
        y0.append(flavor_sign[j]*nu_spec)

    #t array
    t_i = t_i*c_const*from_eV_to_1_over_m #from s to eV⁻¹
    t_f = t_f*c_const*from_eV_to_1_over_m #from s to eV⁻¹
    t_step=(2*np.pi/max(omega,mu_0))/30
    t = np.arange(t_i,t_f,t_step) #eV⁻¹

    return y0,omega,E_nu,t,mu_0,n_f,n_dim,n_theta,cos_theta_vec


def func_FFC_uniform(time,y, params):
    omega,E,theta_V,mu_0,n_f,n_dim,n_theta,cos_theta_vec= params  # unpack parameters
    B=np.array(B_vec(n_dim,theta_V))
    L=np.array(L_vec(n_dim))
    mu=mu_0
    lamb=lambda_supernova(time,"no",0)

    derivs=[]
    nu, nubar = from_1D_to_MultiD(y,n_f,n_dim,n_theta)
    
    #Summed nu and nubar components
    nu_sum, nubar_sum=[],[]
    nu_sum_theta, nubar_sum_theta=[],[]
    nu_aux=np.transpose(nu,(2,0,1))
    nubar_aux=np.transpose(nubar,(2,0,1))

    for i in range(n_dim):
      nu_sum.append(sum(map(sum,nu_aux[i])))
      nubar_sum.append(sum(map(sum,nubar_aux[i])))
      nu_sum_theta.append(sum(cos_theta_vec*list(map(sum,nu_aux[i]))))
      nubar_sum_theta.append(sum(cos_theta_vec*list(map(sum,nubar_aux[i]))))
    
    nu_sum=np.array(nu_sum)
    nubar_sum=np.array(nubar_sum)
    nu_sum_theta=np.array(nu_sum_theta)
    nubar_sum_theta=np.array(nubar_sum_theta)
    
    # list of dy/dt=f functions
    for i in range(n_theta):
      for j in range(n_f):
        #nu
        aux=B*omega+L*lamb-mu*((nu_sum-nubar_sum)-cos_theta_vec[i]*(nu_sum_theta-nubar_sum_theta))
        P_aux= cross_prod(nu[i][j],aux)
        for k in range(n_dim):
          derivs.append(P_aux[k])
        
        #nubar
        aux=-1*B*omega+L*lamb-mu*((nu_sum-nubar_sum)-cos_theta_vec[i]*(nu_sum_theta-nubar_sum_theta))
        P_aux= cross_prod(nubar[i][j],aux)
        for k in range(n_dim):
          derivs.append(P_aux[k])

    return derivs

def solver_FFC_uniform(nu_types,t_i,t_f,E_nu,delta_m2,theta_V,mu_0,theta_bins,case,mass_ord,solver,seed=0,atol=8e-10,rtol=8e-13):
    y0,omega,E,t,mu_0,n_f,n_dim,n_theta,cos_theta_vec=initiate_FFC_uniform(nu_types,t_i,t_f,E_nu,delta_m2,mu_0,theta_bins,case,seed)
    
    if mass_ord=="NH": 
        params=omega,E,theta_V,mu_0,n_f,n_dim,n_theta,np.array(cos_theta_vec)
    elif mass_ord=="IH":
        params=-1*omega,E,theta_V,mu_0,n_f,n_dim,n_theta,np.array(cos_theta_vec)
    else:
        print("Not a mass ordering option!")
    if solver=="odeint":
        y_sol= odeint(func_FFC_uniform, y0, t, args=(params,),tfirst=True,rtol=rtol,atol=atol)
    else:
        psoln= solve_ivp(func_FFC_uniform,t_span=[t[0],t[-1]] ,y0=y0, method=solver, args=(params,),t_eval=t,rtol=rtol,atol=atol)
        y_sol=np.transpose(psoln.y)
        t=psoln.t
    nu, nubar= read_output(y_sol,(n_f,n_dim,n_theta))    
    t=t/(c_const*from_eV_to_1_over_m) #From eV⁻¹ to km

    return cos_theta_vec,t,mu_0, nu, nubar
