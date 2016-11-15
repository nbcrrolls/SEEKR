#!/usr/bin/python

'''
Fixes the VMD output PDB file to have consecutive number and correctly formatted TER cards
'''

import os, sys
import pdb2 as pdb
parser = pdb.Big_PDBParser()


def load_and_save_pdb(infilename, outfilename=''):
  if not outfilename: outfilename = '.'.join(os.path.basename(infilename).split('.')[:-1])+'_fixed.pdb'
  print "outfilename:", outfilename
  pdb = parser.get_structure('pdb', infilename, preserve_resid=False)
  pdb.save(outfilename, amber=True, standard=False)
  print "file saved"
  return
  
if __name__ == "__main__":
  assert len(sys.argv) == 2, "Incorrect number of arguments"
  load_and_save_pdb(sys.argv[1])
