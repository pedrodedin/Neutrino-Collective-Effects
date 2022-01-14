from Auxiliar_Functions import *

def initiate(nu_types,E_i,E_f,E_step,E_0,Amplitude):
	y0=[] #Initial state
	B=[]
	flavor_sign=1

	E_vec=np.arange(E_i,E_f,E_step)
	n_E=len(E_vec)

	n_f=len(nu_types)
	n_dim=(n_f**2)-1

	for i in range(n_E): 
	  B.append([])
	  for j in range(n_dim):     
	    B[i].append(B_vec(E_vec[i],n_dim)[1+j])
	    
	  for j in range(n_f):
	    if nu_types[j]=="nu_x":
	      flavor_sign=-1
	    if nu_types[j]=="nu_e":
	      flavor_sign=1
	    #nu
	    nu_spec=Amplitude[n_f*j]*phi_vec(E_vec[i],E_0[n_f*j+1],2.3)*E_step
	    y0.append(0)
	    y0.append(0)
	    y0.append(flavor_sign*nu_spec)
	    #nubar
	    nu_spec=Amplitude[n_f*j+1]*phi_vec(E_vec[i],E_0[n_f*j+1],2.3)*E_step
	    y0.append(0)
	    y0.append(0)
	    y0.append(flavor_sign*nu_spec)

	#time
	t_bins=1000
	t_max = 4*(2*np.pi/min(B_vec(E_vec,n_dim)[0])) #eV⁻¹
	t_step = (2*np.pi/max(B_vec(E_vec,n_dim)[0]))/20 #eV⁻¹
	t_vec = np.arange(0., t_bins*t_step , t_step) #eV⁻¹

	#mu
	mu_0=(10)*max(B_vec(E_vec,n_dim)[0])

	return y0,B,E_vec,t_vec,mu_0,n_f,n_dim,n_E
