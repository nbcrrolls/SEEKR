import string, re, glob, argparse
import numpy as np
from numpy import linalg as LA
from pprint import pprint
import MDAnalysis as mda
from MDAnalysis import *
from MDAnalysis.analysis.align import *

FORWARD_GLOB= "forward.*.dcd"
#parmfile_name= "../building/holo.parm7"
#equil_filename= '../ens_equil/ens_equil_center.dcd'
#lig_res = "BEN"
#dcd_spacing=1


def collect_forwards():
  traj_number= []
  forward_filenames= glob.glob(FORWARD_GLOB)
  for file in forward_filenames:
    traj_number.append(int(file.split('.')[1]))
  return(traj_number)

def extract_fhpd(lig_res, traj_numbers,parmfile_name, equil_filename, top, dcd_spacing):
  LIG_x = []
  LIG_y = []
  LIG_z = []
  ens_traj= (traj_numbers*dcd_spacing) # using ens_equil_center.dcd which has already removed starting frames
  u = mda.Universe(parmfile_name,equil_filename, topology_format= top)
  for ts in u.trajectory:
    if ts.frame in ens_traj: 
      LIG_x.append(u.select_atoms("resname %s" %lig_res).center_of_mass()[0])
      LIG_y.append(u.select_atoms("resname %s" %lig_res).center_of_mass()[1])
      LIG_z.append(u.select_atoms("resname %s" %lig_res).center_of_mass()[2])
  return(LIG_x, LIG_y, LIG_z)



def write_pdb(LIG_x, LIG_y, LIG_z):
  outfile=open('lig_com_fhpd.pqr','w')

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
  parser = argparse.ArgumentParser(description="A program that extracts the first hitting point distribution from SEEKR forward/reverse calculations.")
  parser.add_argument('struct', metavar='STRUCT', type=str, help='the path to the structure file, for example: ../building/holo.parm7')
  parser.add_argument('-top','--top_format', metavar='TOP_FORMAT', type=str, default='PRMTOP', help='the format of the structure file given above for use by MDAnalysis. Default: PRMTOP')
  parser.add_argument('traj', metavar='TRAJ', type=str, help='the path to the EQUILIBRIUM trajectory file-- should be RMSD aligned')
  parser.add_argument('dcd_space', metavar='DCD_SPACE', type=int, help='the integer value of the dcd spacing used to initialize reversals')
  parser.add_argument('lig', metavar='LIG', type=str, help='the residue name of the ligand molecule in the structure file')

  args = parser.parse_args() # parse all the arguments
  args = vars(args)

  struct= args['struct']
  traj= args['traj']
  lig= args['lig']
  top= args['top_format']
  dcd_spacing= args['dcd_space']


  traj_numbers=collect_forwards()
  print len(traj_numbers), " successful reversals"
  LIG_x, LIG_y, LIG_z =extract_fhpd(lig, traj_numbers, struct, traj, top, dcd_spacing)
  write_pdb(LIG_x, LIG_y, LIG_z)  


if __name__ == "__main__": main()
