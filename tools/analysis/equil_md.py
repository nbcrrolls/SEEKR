import argparse
import string
import numpy as np
from numpy import linalg as LA
from pprint import pprint
import MDAnalysis as mda


def extract_traj_info(struct_filename, traj_center_filename, lig_res, top):
  u = mda.Universe(struct_filename,traj_center_filename,topology_format=top)
  LIG_x = []
  LIG_y = []
  LIG_z = []
  for ts in u.trajectory:
    LIG_x.append(u.select_atoms("resname %s" %lig_res).center_of_mass()[0])
    LIG_y.append(u.select_atoms("resname %s" %lig_res).center_of_mass()[1])
    LIG_z.append(u.select_atoms("resname %s" %lig_res).center_of_mass()[2])

  return(LIG_x, LIG_y, LIG_z)

def write_pdb(LIG_x, LIG_y, LIG_z):
  outfile=open('lig_com_equil.pqr','w')
  
  atom= 'ATOM'
  serial= 1
  atomname= 'CA'
  alt=' ' 
  resname= 'AAA'
  chain= 'B'
  resnumber= 1
  code=' '
  occ= 1.00
  temp = 0.00
  elem = 'C'
  charge= 00
  
  
  for frame in range(len(LIG_x)):
    outfile.write('%-6s%5d %4s%1s%3s %1s%8d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s' % (atom,serial,atomname, alt, resname,chain, resnumber,code,  LIG_x[frame], LIG_y[frame], LIG_z[frame], occ, temp, elem, charge) +'\n') 
    resnumber +=1

def main():
  parser = argparse.ArgumentParser(description="A program that extracts the equilibrium distribution from SEEKR umbrella sampling calculations.")
  parser.add_argument('struct', metavar='STRUCT', type=str, help='the path to the structure file, for example: ../building/holo.parm7')
  parser.add_argument('-top','--top_format', metavar='TOP_FORMAT', type=str, default='PRMTOP', help='the format of the structure file given above for use by MDAnalysis. Default: PRMTOP')
  parser.add_argument('traj', metavar='TRAJ', type=str, help='the path to the trajectory file-- should be RMSD aligned')
  parser.add_argument('lig', metavar='LIG', type=str, help='the residue name of the ligand molecule in the structure file')

  args = parser.parse_args() # parse all the arguments
  args = vars(args)

  struct= args['struct']
  traj= args['traj']
  lig= args['lig']
  top= args['top_format']

 # struct_filename=str('../building/holo.parm7')
 # traj_center_filename=str('ens_equil_center.dcd')
 # lig_res=str('BEN') 
  
  print 'structure file:', struct
  print 'centered trajectory file:', traj
  print 'lig_res', lig 
 
  #print 'structure file:', struct_filename
  #print 'centered trajectory file:', traj_center_filename
  #print 'lig_res', lig_res
  LIG_x, LIG_y, LIG_z = extract_traj_info(struct, traj, lig, top)
  write_pdb(LIG_x, LIG_y, LIG_z)

if __name__ == "__main__": main()
