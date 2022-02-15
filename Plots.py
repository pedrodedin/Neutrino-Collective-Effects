from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
import matplotlib.pyplot as plt
import numpy as np
from Auxiliar_Functions import *

def Plot_Probability_Solar(S3,r_per_R_sol,E):
    Pee_high_E = np.sin(theta_21)**2
    Pee_low_E = 1-(1/2)*np.sin(2*theta_21)**2

    fig = plt.figure(1, figsize=(15,5)) 

    plt.subplot(121)
    plt.plot(r_per_R_sol, 1/2*(1+S3),label=r'E=%.1f MeV'%(E))
    plt.axhline(Pee_high_E,c='r',ls='--',label=r'$\left< \overline{P}_{\nu_e \rightarrow \nu_e}^{High} \right> = sin^2\theta$')
    plt.axhline(Pee_low_E,c='y',ls='--',label=r'$\left< \overline{P}_{\nu_e \rightarrow \nu_e}^{Low} \right> = 1-\frac{1}{2}sin^2 2\theta$')
    plt.xlabel(r'$r/R_{\odot}$')
    plt.ylabel('$P_{ee}$')
    plt.title(r"Survival Probability - P($\nu_e \rightarrow \nu_e$)")
    plt.legend(loc='upper right')

    plt.subplot(122)
    plt.plot(r_per_R_sol, 1/2*(1-S3),label=r'E=%.1f MeV'%(E))
    plt.axhline(1-Pee_high_E,c='r',ls='--',label=r'$\left< \overline{P}_{\nu_e \rightarrow \nu_\mu}^{High} \right> = 1-sin^2\theta$')
    plt.axhline(1-Pee_low_E,c='y',ls='--',label=r'$\left< \overline{P}_{\nu_e \rightarrow \nu_e}^{Low} \right> = \frac{1}{2}sin^2 2\theta$')
    plt.xlabel(r'$r/R_{\odot}$')
    plt.ylabel('$P_{e\mu}$')
    plt.title(r"Conversion Probability - P($\nu_e \rightarrow \nu_\mu$)")
    plt.legend(loc='lower right')

    plt.tight_layout()
    plt.show()

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
    ax2.set_ylabel("$P_{e}$")
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

def Plot_Spectrum(E_vec,E_0,nu_e,nubar_e, nu_x,nubar_x,title):
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
    ax1.set_ylabel(r'$\phi(E) [MeV^{-1}]$')
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
    ax2.set_ylabel(r'$\phi(E) [MeV^{-1}]$')
    ax2.legend(loc='upper right')
    ax2.set_ylim(0,1.1*max(nu_e[0]))
    
    fig.suptitle((r'Isotropic Neutrino Gas - Mass Hierarchy: %s'%(title))+"\n"+
                 (r"$\overline{E}_{\nu_e}= %.1f$ MeV, $\overline{E}_{\overline{\nu}_e}= %.1f$ MeV, $\overline{E}_{\nu_x}= %.1f$ MeV"%(E_0[0],E_0[1],E_0[2])))
    fig.savefig("Figures/Isotrpic_Gas_Spectrum_%s.png"%title)