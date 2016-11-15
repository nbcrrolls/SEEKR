import argparse
import string, re, glob
import numpy as np
from numpy import linalg as LA
from pprint import pprint
import MDAnalysis as mda

LIG_PQR_GLOB= 'lig*.pqr'
lig_res= 'GHO'

def extract_lig_crd():
  LIG_x = []
  LIG_y = []
  LIG_z = []
  pqr_filenames= glob.glob(LIG_PQR_GLOB)
  for file in pqr_filenames:
    u= mda.Universe(file)
    LIG= np.array(u.select_atoms('resname %s' %lig_res).positions)
    LIG_x.append(LIG[0][0])
    LIG_y.append(LIG[0][1])
    LIG_z.append(LIG[0][2])
  return(LIG_x, LIG_y, LIG_z)

def write_pdb(LIG_x, LIG_y, LIG_z):
  outfile=open('bd_lig_com_fhpd.pqr','w')

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
  LIG_x, LIG_y, LIG_z = extract_lig_crd()
  write_pdb(LIG_x, LIG_y, LIG_z)    

if __name__ == "__main__": main()
