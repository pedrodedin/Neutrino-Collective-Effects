from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np

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
    return ani
