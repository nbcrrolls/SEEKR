#!/usr/bin/python
'''
measure_benz_angle.py
Benjamin Jagger
Amaro Lab 2016

Calculate the angle between the ligand and a vector drawn from the COM of the receptor 

'''
import argparse
import string
import numpy as np
from numpy import linalg as LA
from pprint import pprint
import MDAnalysis as mda

np.set_printoptions(threshold='nan')

def extract_traj_info(struct_filename,traj_filename, top, lig_sel_1, lig_sel_2, rec_com_sel):
  u = mda.Universe(struct_filename,traj_filename,topology_format= top)
  angles= []
  for ts in u.trajectory:
    L1=  u.select_atoms(lig_sel_1).center_of_mass()
    L2= u.select_atoms(lig_sel_2).center_of_mass()
    rec_com = u.select_atoms(rec_com_sel).center_of_mass()
    V1= L1-rec_com
    V2= L2-L1
    angle=np.array(np.arccos(np.dot(V2, V1)/(LA.norm(V2)*LA.norm(V1))))
    angle= np.degrees(angle)
    angles.append(angle)
  return angles

def count_orientation(angle):
  face_1= []
  face_1_counts = 0
  face_2= []
  face_2_counts = 0
  plane = []
  plane_counts =0
  
  for angle in theta:
    if angle <90.0:
      face_1.append(angle)
      face_1_counts += 1
    elif angle > 90.0:
      face_2.append(angle)
      face_2_counts +=1
    else:
      plane.append(angle)
      plane_counts +=1
  total_counts = face_1_counts + face_2_counts
  face_1_prob = face_1_counts/float(total_counts)
  face_2_prob = face_2_counts/float(total_counts)
  print 'Primary Face Counts: ' , face_1_counts, 'Probability: ', face_1_prob 
  print 'Secondary Face Counts: ' , face_2_counts, 'Probability: ', face_2_prob
  print 'Total Count: ', total_counts
  print 'Number exactly in plane: ' , plane_counts
  
def write_output(angles):
  outfile= open('lig_orientation.txt', 'w')
  for i in range(len(angles)):
    outfile.write('%d   %d' % (i, angles[i])+ '\n')
  outfile.close()
  outfile2= open('lig_orientation_angles.npy', 'w')
  np.save(outfile2, angles)


def main():
  parser= argparse.ArgumentParser(description='Calculate the angle between a normal vector of a plane of atoms and the COM of ligand')
  parser.add_argument('struct', metavar='STRUCT', type=str, help='the path to the structure file, for example: ../building/holo.parm7')
  parser.add_argument('-top','--top_format', metavar='TOP_FORMAT', type=str, default='PRMTOP', help='the format of the structure file given above for use by MDAnalysis. Default: PRMTOP')
  parser.add_argument('traj', metavar='TRAJ', type=str, help='the path to the trajectory file-- should be RMSD aligned')
  #parser.add_argument('lig', metavar='LIG', type=str, help='the residue name of the ligand molecule in the structure file')
  parser.add_argument('lig_sel_1', metavar='SEL_1', type=str, help='The first atom selection used for determining the ligand vector')
  parser.add_argument('lig_sel_2', metavar='SEL_2', type=str, help='The seconf atom selection used for determining the ligand vector')
  parser.add_argument('rec_com_sel', metavar='REC_COM', type=str, help='The atom selection used for determining the receptor center of mass')


  args=parser.parse_args()
  args=vars(args)
  struct_filename= args['struct']
  traj_filename= args['traj']
  #lig= args['lig']
  top= args['top_format']
  lig_sel_1= args['lig_sel_1']
  lig_sel_2= args['lig_sel_2']
  rec_com_sel= args['rec_com_sel']
  


  print 'structure file:', struct_filename
  print 'trajectory file:', traj_filename 
  print 'first ligand selection:' , lig_sel_1
  print 'second ligand selection:', lig_sel_2
  print 'receptor COM selection:', rec_com_sel
  
  angles =  extract_traj_info(struct_filename,traj_filename, top, lig_sel_1, lig_sel_2, rec_com_sel)
  #pprint(angles)
  write_output(angles)
  #count_orientation(theta)


if __name__ == "__main__": main()
