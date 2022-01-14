from Auxiliar_Functions import *
from Initial_Conditions import *
from scipy.integrate import odeint

def func_Collective_nu(y, time, params):
    B,mu_0,n_f,n_dim,n_E= params  # unpack parameters
    mu=mu_supernova_vec(time,mu_0)
    derivs=[]
    nu, nubar = [],[]
    num_diff_nu_compnents=2*n_f*n_dim

    #Filling [Energy bin][Nu_types][3components]
    for i in range(n_E):
      nu.append([])
      nubar.append([])
      for j in range(n_f):
        nu[i].append([])
        nubar[i].append([])
        for k in range(n_dim):
          #nu 
          nu_index=(i*num_diff_nu_compnents)+k+2*j*n_dim
          nu[i][j].append(y[nu_index])
          #nubar   
          nubar_index=(i*num_diff_nu_compnents)+(k+n_dim)+2*j*n_dim
          nubar[i][j].append(y[nubar_index])
    
    #Summed nu and nubar components
    nu_sum, nubar_sum=[],[]
    nu_aux=np.transpose(nu,(2,0,1))
    nubar_aux=np.transpose(nubar,(2,0,1))

    for i in range(n_dim):
      nu_sum.append(sum(map(sum,nu_aux[i])))
      nubar_sum.append(sum(map(sum,nubar_aux[i])))
    
    B=np.array(B)
    nu_sum=np.array(nu_sum)
    nubar_sum=np.array(nubar_sum)
    
    # list of dy/dt=f functions
    for i in range(n_E):
      for j in range(n_f):
        #nu 
        P_aux= cross_prod(nu[i][j],(B[i]-mu*(nu_sum-nubar_sum)))
        for k in range(n_dim):
          derivs.append(P_aux[k])
        
        #nubar
        P_aux= cross_prod(nu[i][j],(-1*B[i]-mu*(nu_sum-nubar_sum)))
        for k in range(n_dim):
          derivs.append(P_aux[k])

    return derivs


def solver_two_families(nu_types,E_i,E_f,E_step,E_0,Amplitude,mass_ord):

	#E_vec=np.arange(E_i,E_f,E_step)
	y0,B,E_vec,t_vec,mu_0,n_f,n_dim,n_E=initiate(nu_types,E_i,E_f,E_step,E_0,Amplitude)

	if mass_ord=="NH": 
		params=np.array(B),mu_0,n_f,n_dim,n_E
	elif mass_ord=="IH":
		params=-1*np.array(B),mu_0,n_f,n_dim,n_E
	else:
		print("Not a mass ordering option!")
		return 0

	psoln= odeint(func_Collective_nu, y0, t_vec, args=(params,))
	
	nu, nubar= read_output(psoln,params)
	nu_e_time,nubar_e_time,nu_x_time,nubar_x_time=read_two_flavor(nu, nubar)

	#return nu_e_time,nubar_e_time, nu_x_time,nubar_x_time
	return E_vec,t_vec,nu_e_time,nubar_e_time, nu_x_time,nubar_x_time
