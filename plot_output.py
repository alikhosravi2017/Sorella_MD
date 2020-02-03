#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('ggplot')

try:
    with open("output.dat","r") as f:
    	#step,E_kin,E_pot,E_tot,temp_now
    	data = np.loadtxt(f)
    	fig, axs = plt.subplots(2, 1)

    	axs[0].plot(data[:,0],data[:,1],data[:,0],data[:,2],data[:,0],data[:,3])
    	axs[0].set_ylim(top=4*np.min(data[:,1]),bottom=4*np.min(data[:,2]))
    	axs[0].set_xlabel("step [a.u.]")
    	axs[0].set_ylabel("Energy [a.u.]")
    	axs[0].legend([r"$E_\mathrm{kin}$",r"$E_\mathrm{pot}$",r"$E_\mathrm{tot}$"])

    	axs[1].plot(data[:,0],data[:,4])
    	axs[1].set_ylim(bottom=0.)
    	axs[1].set_xlabel("step [a.u.]")
    	axs[1].set_ylabel("Temperature [ K ]")
    	# print(np.std(data[10:,4])) # std deviation of temperature
    	fig.tight_layout()
    	plt.show()
except IOError:
    print("File not accessible")


# doesn't work due to 
# try:
# 	import pytraj as pt
# 	import nglview as nv
# 	traj = pt.load("traj.xyz")
# 	nv.show_pytraj(traj)
# except ModuleNotFoundError:
# 	print("Install NGLView and pytraj")