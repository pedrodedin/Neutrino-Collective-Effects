using SpecialFunctions

############# Constants ########################################
G_F=1.1663787*10^(-23) #eV⁻² - Fermi Constant
delta_m2_31=2.5*10^(-3) #eV² - \Delta m²_31
delta_m2_21=7.53*10^(-5); #eV²
theta_31=asin(sqrt(2.18*10^-2)) #\theta_31
theta_21=asin(sqrt(0.307));
N_A=6.02*10^(23) #Avogadro constant

h_slash=6.582119569*10^-16# eV⋅s
c_const=2.99792458*10^8# m/s
from_eV_to_1_over_m=1/(h_slash*c_const)
from_eV_to_1_over_km=from_eV_to_1_over_m*10^(3)
from_1_over_cm3_to_eV3=(1.973*10^(-5))^3
erg_to_MeV=6.2415*10^5;

############### Matter Effects (2 Families) ########################
#Calculate the effective mass squared difference in matter
# L vector
function L_vec(n_dim)
    if n_dim==3
        return [0,0,1]
    else
        println("Dimension not defined")
    end
end    

# B vector
function B_vec(n_dim,theta)
    if n_dim==3
        B1=-1*sin(2*theta)
        B2=0
        B3=cos(2*theta)
        return [B1,B2,B3]
    else
        println("Dimension not defined")
    end
end;


############# Supernova Info ############################
function phi(E,E_0,alpha)
  N=((alpha+1)^(alpha+1))/(E_0*gamma(alpha+1))
  R=N*((E/E_0)^alpha)*exp((-1)*(alpha+1)*E/E_0)
  return R
end;

#Neutrino potential
function mu_supernova(r,mu_opt,mu_0=0) # r in km
    if mu_opt=="SN"
        R_0=40#km
        if r<R_0
            return mu_0
        end
        return mu_0*(R_0/r)^4 #eV
    
    elseif mu_opt=="const"
        return mu_0
    
    else
        println("Not a mu option!")
        return 0
    end
end;


#Electron density profile
function lambda_supernova(r,option,ne=N_A,t=1) # r in km
  Y_e=0.5
  m_n=1.67492749804*10^-24 #g
  ne_aux=0
  if option=="no"
      ne_aux=0
  elseif option=="const"
      ne_aux=ne
#   elseif option=="SN"
#       ne_aux=(Y_e/m_n)*SN_density_profile(r,t)
  else
    println("Not a matter option")
    return 0
  end
    
  return sqrt(2)*G_F*from_1_over_cm3_to_eV3*ne_aux #eV
end;

################## Read Output ##########################
function from_1D_to_MultiD_Energy(sol,n_f,n_dim,n_E)
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
                push!(nu[i][j],sol[nu_index,:])
                #nubar   
                nubar_index=((i-1)*num_diff_nu_compnents)+(k+n_dim)+2*(j-1)*n_dim
                push!(nubar[i][j],sol[nubar_index,:])
            end
        end
    end
    return nu, nubar
end;

function from_P_vector_to_Spectrum(nu, nubar)
  nu_e_time,nubar_e_time=[],[]
  nu_x_time,nubar_x_time=[],[]

  for l in 1:length(nu[1][1][1]) #time array length
      push!(nu_e_time,[])
      push!(nubar_e_time,[])
      push!(nu_x_time,[])
      push!(nubar_x_time,[])
    
      for i in 1:length(nu) #Energy array length
        #nu
        Pee=(1/2)*(1+nu[i][1][3][l]/nu[i][1][3][1])
        if isnan(Pee)
                Pee=0
        end
        Pxx=(1/2)*(1+nu[i][2][3][l]/nu[i][2][3][1])
        if isnan(Pxx)
                Pxx=0
        end
        push!(nu_e_time[l],Pee*nu[i][1][3][1]+(1-Pxx)*(-1)*nu[i][2][3][1])
        push!(nu_x_time[l],Pxx*(-1)*nu[i][2][3][1]+(1-Pee)*nu[i][1][3][1])

        #nubar
        Pee=(1/2)*(1+nubar[i][1][3][l]/nubar[i][1][3][1])
        if isnan(Pee)
                Pee=0
        end
        Pxx=(1/2)*(1+nubar[i][2][3][l]/nubar[i][2][3][1])
        if isnan(Pxx)
            Pxx=0
        end

        push!(nubar_e_time[l],Pee*nubar[i][1][3][1]+(1-Pxx)*(-1)*nubar[i][2][3][1])
        push!(nubar_x_time[l],Pxx*(-1)*nubar[i][2][3][1]+(1-Pee)*nubar[i][1][3][1])
    end
  end
  return   nu_e_time,nubar_e_time, nu_x_time,nubar_x_time
end;


################### Group ######################
function structure_constant(n_dim,i,j,k)
    if i==j || i==k || j==k
        return 0
    end
    
    if n_dim==3
        i_next=i%3+1
        if  i_next == j
            return 1
        end
        if i_next == k
            return -1 
        end
    else
        println("Dimension not defined")
    end
end;
                
function cross_prod(A,B)
    n_components=length(A)
    C=[]
    for i in 1:n_components
      sum=0
      for j in 1:n_components
          for k in 1:n_components
            sum=sum+structure_constant(n_components,i,j,k)*A[j]*B[k]
          end
      end
    push!(C,sum)
    end
    return C
end;

################## Arrow Plots ########################
using LinearAlgebra
function arrow3d!(x, y, z,  u, v, w; as=0.1, lc=:black, la=1, lw=0.4, scale=:identity,label=false)
    (as < 0) && (nv0 = -maximum(norm.(eachrow([u v w]))))
    for (x,y,z, u,v,w) in zip(x,y,z, u,v,w)
        nv = sqrt(u^2 + v^2 + w^2)
        v1, v2 = -[u,v,w]/nv, nullspace(adjoint([u,v,w]))[:,1]
        v4 = (3*v1 + v2)/3.1623  # sqrt(10) to get unit vector
        v5 = v4 - 2*(v4'*v2)*v2
        (as < 0) && (nv = nv0) 
        v4, v5 = -as*nv*v4, -as*nv*v5
        plot!([x,x+u], [y,y+v], [z,z+w], lc=lc, la=la, lw=lw, scale=scale, label=label)
        plot!([x+u,x+u-v5[1]], [y+v,y+v-v5[2]], [z+w,z+w-v5[3]], lc=lc, la=la, lw=lw, label=false)
        plot!([x+u,x+u-v4[1]], [y+v,y+v-v4[2]], [z+w,z+w-v4[3]], lc=lc, la=la, lw=lw, label=false)
    end
end

function arrow3d(x, y, z,  u, v, w; as=0.1, lc=:black, la=1, lw=0.4, scale=:identity,label=false)
    (as < 0) && (nv0 = -maximum(norm.(eachrow([u v w]))))
    for (x,y,z, u,v,w) in zip(x,y,z, u,v,w)
        nv = sqrt(u^2 + v^2 + w^2)
        v1, v2 = -[u,v,w]/nv, nullspace(adjoint([u,v,w]))[:,1]
        v4 = (3*v1 + v2)/3.1623  # sqrt(10) to get unit vector
        v5 = v4 - 2*(v4'*v2)*v2
        (as < 0) && (nv = nv0) 
        v4, v5 = -as*nv*v4, -as*nv*v5
        plot([x,x+u], [y,y+v], [z,z+w], lc=lc, la=la, lw=lw, scale=scale, label=label)
        plot!([x+u,x+u-v5[1]], [y+v,y+v-v5[2]], [z+w,z+w-v5[3]], lc=lc, la=la, lw=lw, label=false)
        plot!([x+u,x+u-v4[1]], [y+v,y+v-v4[2]], [z+w,z+w-v4[3]], lc=lc, la=la, lw=lw, label=false)
    end
end
