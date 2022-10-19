from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
import matplotlib.pyplot as plt
import numpy as np
from Auxiliar_Functions import *

import os
os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin/'
import matplotlib.colors as mcolors
plt.style.use(['science',"high-vis"])
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=mcolors.TABLEAU_COLORS)
# plt.style.use('default')
# plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 15
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['axes.labelsize'] = 18
plt.rcParams["legend.frameon"]=True

plt.rcParams['lines.linewidth'] = 1.5
# plt.rcParams['xtick.major.size'] = 5.0
# plt.rcParams['xtick.minor.size'] = 3.0
# plt.rcParams['ytick.major.size'] = 5.0
# plt.rcParams['ytick.minor.size'] = 3.0


def Plot_Probability_Solar(S3,r_per_R_sol,E):
    Pee_high_E = np.sin(theta_21)**2
    Pee_low_E = 1-(1/2)*np.sin(2*theta_21)**2

    fig = plt.figure(1, figsize=(15,5)) 

    plt.subplot(121)
    plt.plot(r_per_R_sol, 1/2*(1+S3),label=r'E=%.1f MeV'%(E))
    if E>2:
        plt.axhline(Pee_high_E,c='r',ls='--',label=r'$\left< \overline{P}_{\nu_e \rightarrow \nu_e}^{High} \right> = sin^2\theta$')
    else:
        plt.axhline(Pee_low_E,c='y',ls='--',label=r'$\left< \overline{P}_{\nu_e \rightarrow \nu_e}^{Low} \right> = 1-\frac{1}{2}sin^2 2\theta$')
    plt.xlabel(r'$r/R_{\odot}$')
    plt.ylabel('$P_{ee}$')
    plt.title(r"Survival Probability - P($\nu_e \rightarrow \nu_e$)")
    plt.legend(loc='upper right')

    plt.subplot(122)
    plt.plot(r_per_R_sol, 1/2*(1-S3),label=r'E=%.1f MeV'%(E))
    if E>2:   
        plt.axhline(1-Pee_high_E,c='r',ls='--',label=r'$\left< \overline{P}_{\nu_e \rightarrow \nu_\mu}^{High} \right> = 1-sin^2\theta$')
    else:
        plt.axhline(1-Pee_low_E,c='y',ls='--',label=r'$\left< \overline{P}_{\nu_e \rightarrow \nu_e}^{Low} \right> = \frac{1}{2}sin^2 2\theta$')
    plt.xlabel(r'$r/R_{\odot}$')
    plt.ylabel('$P_{e\mu}$')
    plt.title(r"Conversion Probability - P($\nu_e \rightarrow \nu_\mu$)")
    plt.legend(loc='lower right')

    plt.tight_layout()
    return fig
    #plt.show()

def Pol_Vec_Anim_Solar(S,B,r_per_R_sol,E):
    r_f=len(r_per_R_sol)
    r_step=int(r_f/100)
    if r_step<1:
        r_step=1

    fig = plt.figure(figsize=(16, 5), dpi= 80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(1,2,1,projection='3d')

    B1,B2,B3=B[0][0],B[1][0],B[2][0]
    quiver_B = ax1.quiver(0, 0, 0, B1, B2, B3, arrow_length_ratio=0.05,color='r',label=r'$\vec{B}$',normalize=True)

    P1,P2,P3= S[0][0],S[1][0],S[2][0]
    quiver = ax1.quiver(0, 0, 0,P1,P2,P3, arrow_length_ratio=0.05,color='b',label=r'$\vec{P}$')

    ax1.set_xlim(-1, 1)
    ax1.set_ylim(-1, 1)
    ax1.set_zlim(-1, 1)
    ax1.quiver(-0.5, 0, 0, 1, 0, 0, arrow_length_ratio=0.05,color='k')
    ax1.quiver(0, -0.5, 0, 0, 1, 0, arrow_length_ratio=0.05,color='k')
    ax1.quiver(0, 0, -0.5, 0, 0, 1, arrow_length_ratio=0.05,color='k')
    ax1.set_title(r'Polarization Vectors - $E_\nu=%.1f$ MeV'%(E))
    ax1.legend()

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.set_ylim(0, 1)
    ax2.set_xlim(0, r_per_R_sol[r_f-1])
    ax2.set_xlabel(r'$r/R_{\odot}$')
    ax2.set_ylabel("$P_{ee}$")
    ax2.set_title(r"P($\nu_e \rightarrow \nu_e$)")
    Pe, = ax2.plot([], [])

    def update(t_i):
        #3D plot
        ax1.cla()
        #Bvector
        B1,B2,B3 =B[0][t_i],B[1][t_i],B[2][t_i]
        quiver_B = ax1.quiver(0, 0, 0, B1, B2, B3, arrow_length_ratio=0.05,color='r',normalize=True,label=r'$\vec{B}$')
        #P vector
        P1,P2,P3= S[0][t_i],S[1][t_i],S[2][t_i]
        quiver = ax1.quiver(0, 0, 0,P1,P2,P3, arrow_length_ratio=0.05,color='b',label=r'$\vec{P}$')
        ax1.set_xlim(-1, 1)
        ax1.set_ylim(-1, 1)
        ax1.set_zlim(-1, 1)
        ax1.quiver(-0.5, 0, 0, 1, 0, 0, arrow_length_ratio=0.05,color='k')
        ax1.quiver(0, -0.5, 0, 0, 1, 0, arrow_length_ratio=0.05,color='k')
        ax1.quiver(0, 0, -0.5, 0, 0, 1, arrow_length_ratio=0.05,color='k')
        ax1.set_title(r'Polarization Vectors - $E_\nu=%.1f$ MeV'%(E))
        ax1.legend()
        #Probability plots
        Pe.set_data(r_per_R_sol[0:t_i], 1/2*(1+S[2][0:t_i]))

    #ani = FuncAnimation(fig, update,fargs=(E), frames=np.arange(0,r_f,r_step), interval=100)
    ani = FuncAnimation(fig, update, frames=np.arange(0,r_f,r_step), interval=100)
    plt.close()
    return ani

    
######################### Collective Effects #######################################



################ Isotropic Mono-eneretic #######################

def Plot_Probability_Isotropic(P_nu_NH,P_nu_IH,r,E,omega,mu_0):
    fig = plt.figure(1, figsize=(8,5)) 
    plt.plot(r, 1/2*(1+P_nu_NH[2]/P_nu_NH[2][0]),label='Normal Hierarchy')
    plt.plot(r, 1/2*(1+P_nu_IH[2]/P_nu_NH[2][0]),ls='--',label='Inverted Hierarchy')
    plt.ylabel(r'$P_{\nu_e},P_{\overline{\nu}_e}$')
    plt.title(r'Survival Probability - $E_{\nu}$=%1.f MeV, $\omega$=%1.e eV, $\mu_0$=%1.e eV'%(E,omega,mu_0))
    plt.xlabel('r [km]')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()

def Plot_Probability_Isotropic_2(P_nu,names,r,E,omega,mu_0):
    fig = plt.figure(1, figsize=(8,5)) 
    for i in range(len(P_nu)):
        plt.plot(r, 1/2*(1+P_nu[i][2]/P_nu[i][2][0]),label=names[i][0],ls=names[i][1],lw=names[i][2])
    plt.ylabel(r'$P_{\nu_e},P_{\overline{\nu}_e}$')
    plt.title(r'Survival Probability - $E_{\nu}$=%1.f MeV, $\omega$=%1.e eV, $\mu_0$=%1.e eV'%(E,omega,mu_0))
    plt.xlabel('r [km]')
    plt.ylim(0,1.05)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()
    

def Pol_Vec_Anim_Isotropic(P_nu,P_nubar,B,r,E,mass_ord):
    r_f=len(r)
    r_step=int(r_f/100)
    if r_step<1:
        r_step=1

    fig = plt.figure(figsize=(16, 5), dpi= 80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(1,2,1,projection='3d')

    B1,B2,B3=B[0][0],B[1][0],B[2][0]
    quiver_B = ax1.quiver(0, 0, 0, B1, B2, B3, arrow_length_ratio=0.05,color='r',label=r'$\vec{B}$',normalize=True)

    P1,P2,P3= P_nu[0][0],P_nu[1][0],P_nu[2][0]
    quiver = ax1.quiver(0, 0, 0,P1,P2,P3, arrow_length_ratio=0.05,color='b',label=r'$\vec{P}_\nu$')
    
    P1,P2,P3= P_nubar[0][0],P_nubar[1][0],P_nubar[2][0]
    quiver = ax1.quiver(0, 0, 0,P1,P2,P3, arrow_length_ratio=0.05,color='g',label=r'$\vec{P}_{\overline\nu}$')

    ax1.set_xlim(-1, 1)
    ax1.set_ylim(-1, 1)
    ax1.set_zlim(-1, 1)
    ax1.quiver(-0.5, 0, 0, 1, 0, 0, arrow_length_ratio=0.05,color='k')
    ax1.quiver(0, -0.5, 0, 0, 1, 0, arrow_length_ratio=0.05,color='k')
    ax1.quiver(0, 0, -0.5, 0, 0, 1, arrow_length_ratio=0.05,color='k')
    ax1.set_title(r'Polarization Vectors - '+mass_ord+' - '+str(E)+' MeV')
    ax1.legend()

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.set_ylim(0, 1)
    ax2.set_xlim(0, r[r_f-1])
    ax2.set_xlabel(r'$r [km]$')
    ax2.set_ylabel("$P_{ee}$")
    ax2.set_title(r"P($\nu_e \rightarrow \nu_e$)")
    Pee_nu, = ax2.plot([], [],label=r"$\nu$")
    Pee_nubar, = ax2.plot([], [],label=r"$\overline{\nu}$")
    ax2.legend(loc='upper right')

    def update(t_i):
        #3D plot
        ax1.cla()
        #Bvector
        B1,B2,B3 =B[0][t_i],B[1][t_i],B[2][t_i]
        quiver_B = ax1.quiver(0, 0, 0, B1, B2, B3, arrow_length_ratio=0.05,color='r',normalize=True,label=r'$\vec{B}$')
        #P vector
        P1,P2,P3= P_nu[0][t_i],P_nu[1][t_i],P_nu[2][t_i]
        quiver = ax1.quiver(0, 0, 0,P1,P2,P3, arrow_length_ratio=0.05,color='b',label=r'$\vec{P}_\nu$')
        
        P1,P2,P3= P_nubar[0][t_i],P_nubar[1][t_i],P_nubar[2][t_i]
        quiver = ax1.quiver(0, 0, 0,P1,P2,P3, arrow_length_ratio=0.05,color='g',label=r'$\vec{P}_{\overline{\nu}}$')
        
        ax1.set_xlim(-1, 1)
        ax1.set_ylim(-1, 1)
        ax1.set_zlim(-1, 1)
        ax1.quiver(-0.5, 0, 0, 1, 0, 0, arrow_length_ratio=0.05,color='k')
        ax1.quiver(0, -0.5, 0, 0, 1, 0, arrow_length_ratio=0.05,color='k')
        ax1.quiver(0, 0, -0.5, 0, 0, 1, arrow_length_ratio=0.05,color='k')
        ax1.set_title(r'Polarization Vectors - '+mass_ord+' - '+str(E)+' MeV')
        ax1.legend()
        #Probability plots
        Pee_nu.set_data(r[0:t_i], 1/2*(1+P_nu[2][0:t_i]))
        Pee_nubar.set_data(r[0:t_i], 1/2*(1+P_nubar[2][0:t_i]))

    #ani = FuncAnimation(fig, update,fargs=(E), frames=np.arange(0,r_f,r_step), interval=100)
    ani = FuncAnimation(fig, update, frames=np.arange(0,r_f,r_step), interval=100)
    plt.close()
    return ani

def Pol_Vec_Anim_Isotropic_mu_Profile(P_nu,P_nubar,B,r,E,mu_0,mass_ord):
    r_f=len(r)
    r_step=int(r_f/100)
    if r_step<1:
        r_step=1

    fig = plt.figure(figsize=(20, 5), dpi= 80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(1,3,1,projection='3d')

    B1,B2,B3=B[0][0],B[1][0],B[2][0]
    quiver_B = ax1.quiver(0, 0, 0, B1, B2, B3, arrow_length_ratio=0.05,color='r',label=r'$\vec{B}$',normalize=True)

    P1,P2,P3= P_nu[0][0],P_nu[1][0],P_nu[2][0]
    quiver = ax1.quiver(0, 0, 0,P1,P2,P3, arrow_length_ratio=0.05,color='b',label=r'$\vec{P}_\nu$')
    
    P1,P2,P3= P_nubar[0][0],P_nubar[1][0],P_nubar[2][0]
    quiver = ax1.quiver(0, 0, 0,P1,P2,P3, arrow_length_ratio=0.05,color='g',label=r'$\vec{P}_{\overline\nu}$')

    ax1.set_xlim(-1, 1)
    ax1.set_ylim(-1, 1)
    ax1.set_zlim(-1, 1)
    ax1.quiver(-0.5, 0, 0, 1, 0, 0, arrow_length_ratio=0.05,color='k')
    ax1.quiver(0, -0.5, 0, 0, 1, 0, arrow_length_ratio=0.05,color='k')
    ax1.quiver(0, 0, -0.5, 0, 0, 1, arrow_length_ratio=0.05,color='k')
    ax1.set_title(r'Polarization Vectors - '+mass_ord+' - '+str(E)+' MeV')
    ax1.legend()

    ax2 = fig.add_subplot(1,3,2)
    ax2.set_ylim(0, 1)
    ax2.set_xlim(0, r[r_f-1])
    ax2.set_xlabel(r'$r [km]$')
    ax2.set_ylabel("$P_{ee}$")
    ax2.set_title(r"P($\nu_e \rightarrow \nu_e$)")
    Pee_nu, = ax2.plot([], [],label=r"$\nu$")
    Pee_nubar, = ax2.plot([], [],label=r"$\overline{\nu}$")
    ax2.legend(loc='upper right')
    
    ax3 = fig.add_subplot(1, 3, 3)
    mu=[]
    for r_i in r:
        mu.append(mu_supernova(r_i,"SN",mu_0))
    ax3.plot(r,mu)
    ax3.axhline(y=delta_m2_31/(2*E*10**6), color="red", linestyle="--", label=r'$\omega$ for $E_{\nu}$=%1.f MeV'%(E))
    ax3.set_xlabel(r'$r$ [km]')
    ax3.set_ylabel(r'$\mu$ [eV]')
    ax3.set_yscale('log')
    ax3.set_title(r"$\mu$ Profile Model for Supernova - $\mu_0=$%.2e eV"%(mu_0))
    mu_point, = ax3.plot(r[0],mu_supernova(r[0],"SN",mu_0),'r.',markersize=12,label=r'$\mu(r)$')
    ax3.legend()

    def update(t_i):
        #3D plot
        ax1.cla()
        #Bvector
        B1,B2,B3 =B[0][t_i],B[1][t_i],B[2][t_i]
        quiver_B = ax1.quiver(0, 0, 0, B1, B2, B3, arrow_length_ratio=0.05,color='r',normalize=True,label=r'$\vec{B}$')
        #P vector
        P1,P2,P3= P_nu[0][t_i],P_nu[1][t_i],P_nu[2][t_i]
        quiver = ax1.quiver(0, 0, 0,P1,P2,P3, arrow_length_ratio=0.05,color='b',label=r'$\vec{P}_{\nu}$')
        
        P1,P2,P3= P_nubar[0][t_i],P_nubar[1][t_i],P_nubar[2][t_i]
        quiver = ax1.quiver(0, 0, 0,P1,P2,P3, arrow_length_ratio=0.05,color='g',label=r'$\vec{P}_{\overline{\nu}}$')
        
        ax1.set_xlim(-1, 1)
        ax1.set_ylim(-1, 1)
        ax1.set_zlim(-1, 1)
        ax1.quiver(-0.5, 0, 0, 1, 0, 0, arrow_length_ratio=0.05,color='k')
        ax1.quiver(0, -0.5, 0, 0, 1, 0, arrow_length_ratio=0.05,color='k')
        ax1.quiver(0, 0, -0.5, 0, 0, 1, arrow_length_ratio=0.05,color='k')
        ax1.set_title(r'Polarization Vectors - '+mass_ord+' - '+str(E)+' MeV')
        ax1.legend()
        
        #Probability plots
        Pee_nu.set_data(r[0:t_i], 1/2*(1+P_nu[2][0:t_i]))
        Pee_nubar.set_data(r[0:t_i], 1/2*(1+P_nubar[2][0:t_i]))
        
        #mu profile
        mu_point.set_data(r[t_i],mu_supernova(r[t_i],"SN",mu_0))
        
    #ani = FuncAnimation(fig, update,fargs=(E), frames=np.arange(0,r_f,r_step), interval=100)
    ani = FuncAnimation(fig, update, frames=np.arange(0,r_f,r_step), interval=100)
    plt.close()
    return ani





################ Isotropic Spectrum #######################

def animation_2_families_spectrum(E_vec,t_vec,nu_e,nubar_e, nu_x,nubar_x,title):
    
    t_step=5
    t_f=len(t_vec)

    fig= plt.figure(figsize=(7, 5), dpi= 80, facecolor='w', edgecolor='k')
    ax1=fig.add_subplot(111)
    ax1.set_xlabel(r'$E [MeV]$')
    ax1.set_ylabel(r'$\phi(E) [MeV^{-1}]$')
    #nu_e
    nu_line, = ax1.plot(E_vec, nu_e[0],color='b',label=r'$\nu_e$')
    ax1.plot(E_vec, nu_e[0],color='b', linestyle="--",label=r'$\nu_e(t=0)$')#Initial
    #nubar_e
    antinu_line, = ax1.plot(E_vec,nubar_e[0],color='r',label=r'$\overline{\nu}_e$')
    ax1.plot(E_vec, nubar_e[0],color='r', linestyle="--",label=r'$\overline{\nu}_e(t=0)$')#Initial
    
    #nu_x
    nu_x_line, = ax1.plot(E_vec,nu_x[0],color='g',label=r'$\nu_x$')
    ax1.plot(E_vec, nu_x[0],color='g', linestyle="--",label=r'$\nu_x(t=0)$')#Initial
    #nubar_x
    antinu_x_line, = ax1.plot(E_vec,nubar_x[0],color='orange',label=r'$\overline{\nu}_x$')
    ax1.plot(E_vec, nubar_x[0],color='orange', linestyle="--",label=r'$\overline{\nu}_x(t=0)$')#Initial
    
    #ax1.set_ylim(0,0.07)
    ax1.set_title(title)
    ax1.legend(loc='upper right')

    def update(t_i):
        nu_line.set_data(E_vec,nu_e[t_i])
        antinu_line.set_data(E_vec,nubar_e[t_i])
        nu_x_line.set_data(E_vec,nu_x[t_i])
        antinu_x_line.set_data(E_vec,nubar_x[t_i])

    ani = FuncAnimation(fig, update, frames=np.arange(0,t_f,t_step), interval=100)
    plt.close()
    return ani

def Plot_Spectrum(E_vec,E_0,mu_0,nu_e,nubar_e, nu_x,nubar_x,title):
    fig= plt.figure(figsize=(15, 6), dpi= 80, facecolor='w', edgecolor='k')
    
    ax1=fig.add_subplot(1,2,1)
    #nu_e
    ax1.plot(E_vec, nu_e[-1],color='b',label=r'$\nu_e$')
    ax1.plot(E_vec, nu_e[0],color='b', linestyle="--",label=r'$\nu_e(t=0)$')#Initial
    #nu_x
    ax1.plot(E_vec,nu_x[-1],color='g',label=r'$\nu_x$')
    ax1.plot(E_vec, nu_x[0],color='g', linestyle="--",label=r'$\nu_x(t=0)$')#Initial
    #Text
    ax1.set_title("Neutrinos")
    ax1.set_xlabel(r'$E [MeV]$')
    ax1.set_ylabel(r'$\phi(E) [a.u.]$')
    ax1.legend(loc='upper right')
    ax1.set_ylim(0,1.1*max(nu_e[0]))
    
    ax2=fig.add_subplot(1,2,2)
    #nubar_e
    ax2.plot(E_vec,nubar_e[-1],color='r',label=r'$\overline{\nu}_e$')
    ax2.plot(E_vec, nubar_e[0],color='r', linestyle="--",label=r'$\overline{\nu}_e(t=0)$')#Initial
    #nubar_x
    ax2.plot(E_vec,nubar_x[-1],color='orange',label=r'$\overline{\nu}_x$')
    ax2.plot(E_vec, nubar_x[0],color='orange', linestyle="--",label=r'$\nu_x(t=0)$')#Initial
    #Text
    ax2.set_title("Antineutrinos")
    ax2.set_xlabel(r'$E [MeV]$')
    ax2.set_ylabel(r'$\phi(E) [a.u.]$')
    ax2.legend(loc='upper right')
    ax2.set_ylim(0,1.1*max(nu_e[0]))
    
    fig.suptitle((r'Isotropic Neutrino Gas - Mass Hierarchy: %s'%(title))+"\n"+
                 (r"$\overline{E}_{\nu_e}= %.1f$ MeV, $\overline{E}_{\overline{\nu}_e}= %.1f$ MeV, $\overline{E}_{\nu_x}= %.1f$ MeV, $\mu_0$=%1.e eV"%(E_0[0],E_0[1],E_0[2],mu_0)))
    plt.tight_layout()
    fig.savefig("Figures/Isotrpic_Gas_Spectrum_%s.png"%title)