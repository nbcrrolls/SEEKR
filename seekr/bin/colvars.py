#!/usr/bin/python

'''
colvars.py

creates the necessary files to implement colvars in NAMD

'''
import sys, os, math, shutil, subprocess #, make_fxd
import unittest
import numpy as np
import pdb2 as pdb
from adv_template import *
import md

verbose = True

# given information needed for a colvar:
#  - location to write files
#  - type of colvar (spherical, planar, etc...)
#  - colvar parameters:
#    > force constant
#    > centers (equilibrium distance)
# will create the colvar file called by NAMD as well as the necessary pdb files

self_path = os.path.dirname(os.path.realpath(__file__)) # get the path to this script

COLVAR_PDB_NAME='colvar.pdb'
COLVAR_SCRIPT_NAME = 'colvar.script'
COLVAR_SPHERE_TEMPLATE = os.path.join(self_path, 'colvar_sphere.template')
COLVAR_PLANEZ_TEMPLATE = os.path.join(self_path, 'colvar_plane.template')

parser = pdb.Big_PDBParser()

def create_colvar_input(directory, colvar_type, params={'TEMPLATE_centers':'0.0', 'TEMPLATE_force':'2.0'}):
  '''creates the colvar input file for use in a namd simulation.
  Input:
  directory: the location to place all the files
  colvar_type: the type of colvar to create:
    - 'sphere': for spherical milestones
    - 'plane': for planar milestones
    - ...
  params: additional parameters needed for colvars, like force, centers, etc. Depends on the colvar_type

  '''
  if colvar_type == "sphere": # then implement the spherical colvars
    use_template = COLVAR_SPHERE_TEMPLATE
  elif colvar_type == "plane": # planar colvars
    use_template = COLVAR_PLANEZ_TEMPLATE

  colvars_script = File_template(use_template,params)
  colvars_script_filename = os.path.join(directory, COLVAR_SCRIPT_NAME)
  colvars_script.save(colvars_script_filename)
  return


def create_colvar_pdb(directory, pdb, groups={():1.00,}):
  ''' creates the colvar pdb file
  Input:
  directory: the location to place all the files
  pdb: pdb file to change occupancy and save with. NOTE: must be deepcopied if you care about the previous occupancies
  groups: a dictionary whose keys are tuples, and whose values are what must be placed into the occupancy column

  '''
  #print "colvar GROUPS:", groups
  for atom in pdb.atoms:
    index = atom.index
    found_index = False
    for key in groups:
      if index in key:
        found_index = True
        atom.occupancy = groups[key]
    if not found_index:
      atom.occupancy = 0.00 # by default

  colvar_pdb_filename = os.path.join(directory, COLVAR_PDB_NAME)
  pdb.save(colvar_pdb_filename, amber=True, standard=False)
  return

def main(settings):
  '''generates the files needed to use colvars in NAMD'''
  i = settings['index']
  colvar_settings = settings['ensemble_equil_settings']['colvar_settings']
  directory = settings['md_file_paths'][i]['ens_equil']
  pdb = settings['configs'][i]
  milestone = settings['milestone_pos_rot_list'][i][0]
  colvar_type = milestone.shape
  if colvar_type == "sphere":
    colvar_center = milestone.dimensions['radius'] # unpack all the necessary values from settings
    params = {}
  elif colvar_type == "plane":
    colvar_center = milestone.dimensions['distance']
    params = {'TEMPLATE_axis':milestone.dimensions['normal']}
  colvar_force = colvar_settings['colvar_force']
  groups={}
  group1_atoms = tuple(md.parse_selection(colvar_settings['colvar_ligand_indeces'])) #tuple(colvar_settings['colvar_ligand_indeces'])
  group2_atoms = tuple(md.parse_selection(colvar_settings['colvar_receptor_indeces'])) #tuple(colvar_settings['colvar_receptor_indeces'])
  groups[group1_atoms] = 1.00 # ligand atoms
  groups[group2_atoms] = 2.00 # receptor atoms
  params.update({'TEMPLATE_centers':colvar_center, 'TEMPLATE_force':colvar_force, 'TEMPLATE_trajfreq': colvar_settings['colvarstrajfrequency'] , 'TEMPLATE_restartfreq': colvar_settings['colvarsrestartfrequency']})

  create_colvar_input(directory, colvar_type, params) # create colvar input file
  create_colvar_pdb(directory, pdb, groups)
  return



class Test_colvars(unittest.TestCase):
  # several test cases to ensure the functions in this module are working properly
  def test_main(self): # test this function
    pass

  def test_create_colvar_input(self):
    test_dir = '/tmp'
    test_type = 'sphere'
    params = {'TEMPLATE_force':'1.9', 'TEMPLATE_centers':'10.0'}
    create_colvar_input(test_dir, test_type, params) # test this function
    self.assertTrue(os.path.exists(os.path.join(test_dir,COLVAR_SCRIPT_NAME))) # verify existence of script
    # test planar colvars

  def test_create_colvar_pdb(self):
    test_dir = '/tmp'
    test_holo = parser.get_structure('small','../test/test_tiny.pdb')
    test_pdb = test_holo
    groups = {(1,2,3):1.0, (4,5,6):2.0}
    create_colvar_pdb(test_dir, test_pdb, groups) # test this function
    self.assertTrue(os.path.exists(os.path.join(test_dir,COLVAR_SCRIPT_NAME))) # verify existence of script



if __name__ == "__main__":
  print "Running unit tests for colvar.py"
  unittest.main()
