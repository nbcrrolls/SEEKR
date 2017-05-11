#!/usr/bin/python

'''
amber_colvars.py

creates the necessary files to implement colvars in AMBER

'''
import sys, os, math, shutil, subprocess #, make_fxd
import unittest
import numpy as np
import pdb2 as pdb
from adv_template import *
import md

verbose = True


# NOTE: Currently amber can only handle distance based restraints/colvars with the pmemd##

# given information needed for a colvar:
#  - location to write files
#  - type of colvar (spherical, planar, etc...)
#  - colvar parameters:
#    > force constant
#    > centers (equilibrium distance)
# will create the amber RST file

self_path = os.path.dirname(os.path.realpath(__file__)) # get the path to this script
RST_template= 'amber_RST.template'
RST_filename= 'COM.RST'
#COLVAR_PDB_NAME='colvar.pdb'
#COLVAR_SCRIPT_NAME = 'colvar.script'
#COLVAR_SPHERE_TEMPLATE = os.path.join(self_path, 'colvar_sphere.template')
#COLVAR_PLANEZ_TEMPLATE = os.path.join(self_path, 'colvar_plane.template')


parser = pdb.Big_PDBParser()

def create_colvar_input(directory, colvar_type, colvar_centers, colvar_force ,group1_atoms, group2_atoms):
  '''creates the RST  for use in a namd simulation.
  Input:
  directory: the location to place all the files
  colvar_type: the type of colvar to create:
    - 'sphere': for spherical milestones
    - 'plane': for planar milestones
    - ...
  '''
  if colvar_type == "sphere": # then implement the spherical colvars
    use_template = RST_template
    print colvar_centers
    r1= str(max(float(colvar_centers)-5.0,0)) #setting the minimium distance for the parabolic potential-- cannot be smaller than 0
    r4=str(float(colvar_centers)+5.0) #setting maximum distance for parabolic potential
    params={'r1':r1, 'r2':colvar_centers, 'r4':r4, 'igr1':str(group1_atoms).strip('[]'), 'igr2':str(group2_atoms).strip('[]'), 'rk2':colvar_force, }
    
  elif colvar_type == "plane": # planar colvars not yet implemented for amber
    raise Exception, "plane restraints not yet implemented for amber package simulations"
    #use_template = RST_template

  RST_file = File_template(use_template,params) #write RST file from template
  colvars_script_filename = os.path.join(directory, RST_filename)
  RST_file.save(colvars_script_filename)
  return


def main(settings):
  '''generates the files needed to use colvars restraints in amber'''
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
  group1_atoms = md.parse_selection(colvar_settings['colvar_ligand_indeces']) #tuple(colvar_settings['colvar_ligand_indeces'])
  group2_atoms = md.parse_selection(colvar_settings['colvar_receptor_indeces']) #tuple(colvar_settings['colvar_receptor_indeces'])
  #groups[group1_atoms] = 1.00 # ligand atoms
  #groups[group2_atoms] = 2.00 # receptor atoms
  #params.update({'TEMPLATE_centers':colvar_center, 'TEMPLATE_force':colvar_force})

  create_colvar_input(directory, colvar_type, colvar_center, colvar_force, group1_atoms, group2_atoms) # create colvar input file
  #create_colvar_pdb(directory, pdb, groups)
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
