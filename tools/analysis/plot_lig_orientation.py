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
import argparse


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
#  print len(angles)
  frame=np.arange(0, len(angles))
  frame=frame*.02
  #print frame
  #print angles
  return(frame, angles)

def main():
  parser= argparse.ArgumentParser(description='plot ligand orientation angle for each milestone')
  parser.add_argument('ligname', type= str, help= 'ligand name for plot title')
  parser.add_argument('data_path', type = str, help='the generic path to find the numpy orientation data eg. "anchor*/md/ens_equil/lig_orientation_angles.npy"') 
  parser.add_argument('x_dim', type= int, help='x dimension of subplots')
  parser.add_argument('y_dim', type= int, help='y dimension of subplots')
  parser.add_argument('save_name', type= str, help='path and name to save the plot')
  args=parser.parse_args()
  args=vars(args)
  ligname= args['ligname']
  x_dim= args['x_dim']
  y_dim= args['y_dim']
  save_name = args['save_name']
  data_path = args['data_path']

  print data_path
  fig, axs =plt.subplots(x_dim,y_dim, sharex='all',sharey='all')
  axs = axs.ravel()
  filenames=glob.glob(data_path)
  filenames_sort=sorted(filenames, key = lambda x: int(x.split('_')[5]))
  print filenames_sort
  for i in range(len(filenames)):
   # filename= 'anchor_'+str(i)+'*/md/ens_equil/benz_orientation_angles.npy'
    filename=filenames_sort[i]
#    print filename
    title_str= filename.split('_')[2].split('/')[1]+' '+ filename.split('_')[3]
#    print title_str
    frame, angle = extract_data(filename)
    axs[i].plot(frame, angle, color='k', linewidth=0.3)
    axs[i].axhline(y=90,color='r', linestyle='--')
#    milestones= [1, 1.5, 2, 3, 4, 6, 8, 10, 12]
#    title_str= str(milestones[i])+ '$\AA$ Milestone'
    axs[i].set_title(title_str)
    axs[i].tick_params(labelsize=8)
    axs[i].set_ylim([0, 200])
    #axs[i].spines["right"].set_visible(false)
  #for ax in axs:  
  #  ax.set_xlabel('Simulation Time (ns)')
  #  ax.set_ylabel('Pocket Volume ($\AA^3$)')
  
  fig.text(0.5, 0.04, 'Simulation Time (ns)', ha='center',)
  fig.text(0.02, 0.5, 'Orientation Angle (degrees)', va='center', rotation='vertical')
  plt.suptitle(ligname + ' Equilibrium Orientations ', fontsize=18)
  plt.tight_layout()
  plt.subplots_adjust(top=0.85, left=0.10, bottom=0.10)
  plt.savefig(save_name, format='png', dpi=500)
#  plt.show()

 
if __name__ == "__main__": main() 
