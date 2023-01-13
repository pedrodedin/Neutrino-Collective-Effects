include("Auxiliar_Functions.jl")
using DifferentialEquations

################ Mono-Energetic ################

function func_Isotropic_Monoenergetic!(dy,y, params, time)
    omega,mu_opt,mu_0,lamb_opt,lamb_0,n_dim= params  # unpack parameters
    
    B=B_vec(n_dim,theta_31)
    L=L_vec(n_dim)
    
    r=time/from_eV_to_1_over_km #From eV⁻¹ to km
    mu=mu_supernova(r,mu_opt,mu_0)
    lamb=lambda_supernova(r,lamb_opt,lamb_0) 
    
    #nu
    P_aux= cross_prod(y[1:n_dim],(B*omega+L*lamb+mu*y[n_dim+1:end]))
    for k in 1:n_dim
        dy[k]=P_aux[k]
    end 
    #nu_bar 
    P_aux= cross_prod(y[n_dim+1:end],(-1*B*omega+L*lamb-mu*y[1:n_dim]))
    for k in 1:n_dim
        dy[k+3]=P_aux[k]
    end 
end;

function solver_Isotropic_Monoenergetic(P,E,r_i,r_f,mass_ord,mu_opt,mu_0;lamb_opt="no",lamb_0=0,n_f=2,reltol=1e-8,abstol=1e-8)

    omega=delta_m2_31/(2*E*10^6) #eV
    r_step = (2*pi/max(omega,mu_0))/200 #eV⁻¹
    r_i = r_i*from_eV_to_1_over_km
    r_f = r_f*from_eV_to_1_over_km
    r = range(r_i,r_f,step=r_step) #eV⁻¹
    n_dim=(n_f^2)-1
    
    if mass_ord=="NH"
        params=omega,mu_opt,mu_0,lamb_opt,lamb_0,n_dim
    elseif mass_ord=="IH"
        omega=-1*omega
        params=omega,mu_opt,mu_0,lamb_opt,lamb_0,n_dim
    else
        println("Not a mass ordering option!")
        return 0
    end

    prob = ODEProblem(func_Isotropic_Monoenergetic!,P,(r_i,r_f),params,reltol=reltol,abstol=abstol)
    sol = solve(prob);  
    
    return sol
            
end;

################ Multi-Energetic ################
function initiate_Isotropic_Multienergetic(nu_types,r_i,r_f,E_i,E_f,E_step,E_0,Amplitude)
    y0=Float64[]#Initial state
    omega=Float64[]
    flavor_sign=1

    E_vec=range(E_i,E_f,step=E_step)
    n_E=length(E_vec)

    n_f=length(nu_types)
    n_dim=(n_f^2)-1
    
    for i in 1:n_E 
        push!(omega,delta_m2_31/(2*E_vec[i]*10^6)) #eV        
        for j in 1:n_f
            if nu_types[j]=="nu_x"
              flavor_sign=-1
            end
            if nu_types[j]=="nu_e"
              flavor_sign=1
            end
        #nu
        index=(n_f*(j-1))+1
        nu_spec=Amplitude[index]*phi(E_vec[i],E_0[index],2.3)*E_step
        push!(y0,0.0)
        push!(y0,0.0)
        push!(y0,flavor_sign*nu_spec)
        #nubar
        nu_spec=Amplitude[index+1]*phi(E_vec[i],E_0[index+1],2.3)*E_step
        push!(y0,0.0)
        push!(y0,0.0)
        push!(y0,flavor_sign*nu_spec)
        end
    end

    #mu
    mu_0=(10)*maximum(omega)
    #r array
    r_step = (2*pi/maximum(omega))/20 #eV⁻¹
    r_i = r_i*from_eV_to_1_over_km #eV⁻¹
    r_f = r_f*from_eV_to_1_over_km #eV⁻¹
    r = range(r_i,r_f,step=r_step) #eV⁻¹

    return y0,omega,E_vec,r,mu_0,n_f,n_dim,n_E

end;

function func_Isotropic_Multienergetic!(dy,y, params, time)
    omega,mu_opt,mu_0,n_f,n_dim,n_E= params  # unpack parameters
    
    B=B_vec(n_dim,theta_31)
    L=L_vec(n_dim)
    
    r=time/from_eV_to_1_over_km #From eV⁻¹ to km
    mu=mu_supernova(r,mu_opt,mu_0)
    lamb=lambda_supernova(r,"no",0) 
    
    nu, nubar = [],[]
    num_diff_nu_compnents=2*n_f*n_dim
    #Filling [Energy bin][Nu_types][3components]
    for i in 1:n_E
        push!(nu,[])
        push!(nubar,[])
        for j in 1:n_f
            push!(nu[i],[])
            push!(nubar[i],[])
            for k in 1:n_dim
                #nu 
                nu_index=((i-1)*num_diff_nu_compnents)+k+2*(j-1)*n_dim
                push!(nu[i][j],y[nu_index])
                #nubar   
                nubar_index=((i-1)*num_diff_nu_compnents)+(k+n_dim)+2*(j-1)*n_dim
                push!(nubar[i][j],y[nubar_index])
            end
        end
    end
    
    #Summed nu and nubar components
    nu_sum=sum(sum(nu))
    nubar_sum=sum(sum(nubar))

    # list of dy/dt=f functions
    for i in 1:n_E
        for j in 1:n_f
            #nu 
            P_aux= cross_prod(nu[i][j],(B*omega[i]+L*lamb-mu*((nu_sum-nu[i][j])-nubar_sum)))
            for k in 1:n_dim
              nu_index=((i-1)*num_diff_nu_compnents)+k+2*(j-1)*n_dim
              dy[nu_index]=P_aux[k]
            end
            #nubar
            P_aux= cross_prod(nubar[i][j],(-1*B*omega[i]+L*lamb-mu*(nu_sum-(nubar_sum-nubar[i][j]))))
            for k in 1:n_dim
              nubar_index=((i-1)*num_diff_nu_compnents)+(k+n_dim)+2*(j-1)*n_dim
              dy[nubar_index]=P_aux[k]
            end
        end
    end
    
end;


function solver_Isotropic_Multienergetic(nu_types,r_i,r_f,E_i,E_f,E_step,E_0,Amplitude,mass_ord;reltol=1e-8,abstol=1e-8)

    y0,omega,E_vec,r,mu_0,n_f,n_dim,n_E=initiate_Isotropic_Multienergetic(nu_types,r_i,r_f,E_i,E_f,E_step,E_0,Amplitude)
    
    if mass_ord=="NH"
        params=omega,"SN",mu_0,n_f,n_dim,n_E
    elseif mass_ord=="IH"
        params=-1*omega,"SN",mu_0,n_f,n_dim,n_E
    else
        print("Not a mass ordering option!")
        return 0
    end
    
    r_i=r_i*from_eV_to_1_over_km
    r_f=r_f*from_eV_to_1_over_km
    
#     func_Collective_nu!(y0,y0,params,1)
    prob = ODEProblem(func_Isotropic_Multienergetic!,y0,(r_i,r_f),params,reltol=reltol,abstol=abstol)
    sol = solve(prob)
    
    return sol,n_f,n_dim,n_E,E_vec,sol.t/from_eV_to_1_over_km

end;
