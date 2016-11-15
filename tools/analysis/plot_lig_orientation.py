#!/usr/bin/python
'''
plot_lig_orientation.py
Benjamin Jagger
Amaro Lab 2016

Plot the angle between a normal vector of a plane of atoms and the COM of ligand-- uses the .npy file from measure_angle.py

'''



import numpy as np
import matplotlib.pyplot as plt
import glob


DATA_GLOB= 'anchor*/md/ens_equil/lig_orientation_angles.npy'

def plot_single(frame, volume):
  plt.plot(frame, volume, 'k')
  plt.ylabel('Pocket Volume ($\AA^3$)')
  plt.xlabel('Simulation Time (ns)')

  plt.show()
  return



def extract_data(filename):
  volume=[]
  angles=np.load(filename)
  #print len(angles)
  frame=np.arange(0, len(angles))
  frame=frame*.02
  #print frame
  #print angles
  return(frame, angles)

def main():
  fig, axs =plt.subplots(2,3, sharex='all',sharey='all')
  axs = axs.ravel()
  filenames=sorted(glob.glob(DATA_GLOB))
  print filenames
  for i in range(len(filenames)):
   # filename= 'anchor_'+str(i)+'*/md/ens_equil/benz_orientation_angles.npy'
    filename=filenames[i]
    print filename
    frame, angle = extract_data(filename)
    axs[i].plot(frame, angle, color='k', linewidth=0.3)
    axs[i].axhline(y=90,color='r', linestyle='--')
    title_str= str((i*2)+2)+ '$\AA$ Milestone'
    axs[i].set_title(title_str)
    axs[i].tick_params(labelsize=8)
    axs[i].set_ylim([0, 200])
    #axs[i].spines["right"].set_visible(false)
  #for ax in axs:  
  #  ax.set_xlabel('Simulation Time (ns)')
  #  ax.set_ylabel('Pocket Volume ($\AA^3$)')
  
  fig.text(0.5, 0.004, 'Simulation Time (ns)', ha='center',)
  fig.text(0.004, 0.5, 'Orientation Angle (degrees)', va='center', rotation='vertical')
  plt.suptitle('Ligand Equilibrium Orientations ', fontsize=18)
  plt.tight_layout()
  plt.subplots_adjust(top=0.85, left=0.10, bottom=0.10)
  plt.savefig('/home/bjagger/images/tryp_paper/lig_angle_plot.png', format='png', dpi=500)
  plt.show()

 
if __name__ == "__main__": main() 
