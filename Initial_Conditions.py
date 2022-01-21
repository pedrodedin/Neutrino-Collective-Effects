from Auxiliar_Functions import *

def initiate(nu_types,t_bins,E_i,E_f,E_step,E_0,Amplitude):
	y0=[] #Initial state
	omega=[]
	flavor_sign=1

	E_vec=np.arange(E_i,E_f,E_step)
	n_E=len(E_vec)

	n_f=len(nu_types)
	n_dim=(n_f**2)-1
    
	for i in range(n_E): 
	  omega.append(delta_m2_31/(2*E_vec[i]*10**6)) #eV        
	  for j in range(n_f):
	    if nu_types[j]=="nu_x":
	      flavor_sign=-1
	    if nu_types[j]=="nu_e":
	      flavor_sign=1
	    #nu
	    nu_spec=Amplitude[n_f*j]*phi_vec(E_vec[i],E_0[n_f*j],2.3)*E_step
	    y0.append(0)
	    y0.append(0)
	    y0.append(flavor_sign*nu_spec)
	    #nubar
	    nu_spec=Amplitude[n_f*j+1]*phi_vec(E_vec[i],E_0[n_f*j+1],2.3)*E_step
	    y0.append(0)
	    y0.append(0)
	    y0.append(flavor_sign*nu_spec)

	#mu
	mu_0=(10)*max(omega)
    
	#time
	t_max = 4*(2*np.pi/min(omega)) #eV⁻¹
	w_max=max(mu_0,max(omega))
	t_step = (2*np.pi/w_max)/5 #eV⁻¹
	t_vec = np.arange(0., t_bins*t_step , t_step) #eV⁻¹

	return y0,omega,E_vec,t_vec,mu_0,n_f,n_dim,n_E

def initiate_v2(nu_types,t_bins,E_i,E_f,E_step,E_0,Amplitude):
	y0=[] #Initial state
	omega=[]
	flavor_sign=1

	E_vec=np.arange(E_i,E_f,E_step)
	n_E=len(E_vec)

	n_f=len(nu_types)
	n_dim=(n_f**2)-1
    
	for i in range(n_E): 
	  omega.append(delta_m2_31/(2*E_vec[i]*10**6)) #eV        
	  #nu
	  nu_e_spec=Amplitude[0]*phi_vec(E_vec[i],E_0[0],2.3)*E_step
	  nu_x_spec=Amplitude[2]*phi_vec(E_vec[i],E_0[2],2.3)*E_step
	  #Pz=(nu_e_spec-nu_x_spec)/(nu_e_spec+nu_x_spec)
	  Pz=(nu_e_spec-nu_x_spec)
	  y0.append(0)
	  y0.append(0)
	  y0.append(Pz)
	  #nubar
	  nu_e_spec=Amplitude[1]*phi_vec(E_vec[i],E_0[1],2.3)*E_step
	  nu_x_spec=Amplitude[3]*phi_vec(E_vec[i],E_0[3],2.3)*E_step
	  #Pz=(nu_e_spec-nu_x_spec)/(nu_e_spec+nu_x_spec)
	  Pz=(nu_e_spec-nu_x_spec)
	  y0.append(0)
	  y0.append(0)
	  y0.append(Pz)
    
	#mu
	mu_0=(10)*max(omega)
    
	#time
	t_max = 4*(2*np.pi/min(omega)) #eV⁻¹
	w_max=max(mu_0,max(omega))
	t_step = (2*np.pi/w_max)/5 #eV⁻¹
	t_vec = np.arange(0., t_bins*t_step , t_step) #eV⁻¹


	return y0,omega,E_vec,t_vec,mu_0,n_f,n_dim,n_E
