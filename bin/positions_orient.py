#!/usr/bin/python

'''
positions_orient.py

Contains functions to be called by seekr.py that generate structural positions/orientations to be merged with a receptor structure.

'''

import pdb2 as pdb # custom library for reading/writing pdb files
import transformations # library for useful rotation functions
import numpy # numerical python
import numpy as np
from numpy import arange
from copy import deepcopy # needed to keep track of separate structure objects
import time, random
from math import sqrt, sin, cos, ceil, floor, log10, atan2, pi, asin
from MDAnalysis import *
from random import random
from subprocess import call
import os
import unittest

#x = numpy.identity(3) # FLAGGED FOR REMOVAL

verbose = True # whether stuff is printed in detail

SO3_grid_program = '/home/lvotapka/Downloads/SO3_grid/SO3_Grid' # NOTE: this should be defined in a configuration file

class Location():
  '''contains the x, y, z coordinates of a ligand center of mass, plus a rotation
matrix to allow further rotation in the case of a polymer receptor'''
  def __init__(self):
    self.x = 0.0
    self.y = 0.0
    self.z = 0.0
    self.matrix = numpy.matrix(numpy.identity(4)) # by default assign the identity rotation to matrix

  def transform(self, matrix):
    coord_vec = numpy.matrix([[self.x], [self.y], [self.z]], 1.0)
    new_vec = matrix * coord_vec # our new rotated coordinates
    self.x = new_vec[0,0]
    self.y = new_vec[1,0]
    self.z = new_vec[2,0]
    self.matrix = matrix * self.matrix

def normalize (quats):
  newquats = []
  for quat in quats:
    newquat = quat / np.linalg.norm(quat)
    newquats.append(newquat.tolist())
  return newquats

def make_SO3_grid(resolution=1):
  hedron = []

  curdir = os.path.realpath('.') # cd into the working directory
  os.chdir('/tmp') # NOTE: this should be customizable
  cmd = ' '.join([SO3_grid_program,'1',str(resolution)])
  call(cmd, shell=True, executable="/bin/bash") # run the SO3 program
  data_file = open('data.qua','r') # read the output
  for line in data_file:
    hedron.append(line.strip().split())
  data_file.close()
  os.chdir(curdir) # cd back to the directory we were at before
  return hedron

def quaternions_to_upper_hemisphere(quat_list, ref_quat=[1,0,0,0]):
  new_quat_list = []
  a_ref_quat = np.array(ref_quat)
  for quat in quat_list:
    a_quat = np.array(quat)
    if np.dot(a_quat,a_ref_quat) >= 0.0:
      new_quat_list.append(quat)
    else:
      new_quat_list.append((-1.0 * a_quat).tolist())
  return new_quat_list

def even_permutation(A):
  '''calculates all the even permutations of an array A
  NOTE: list A must contain unique elements for this algorithm to work
  '''
  maxiter = 500
  numiter = 0
  N = len(A)
  even_perms = [A[:]]
  parityflip = False # keeps track of whether the current permutation is even or odd
  A.sort()
  while True: # iterate thru all even permutations
    while True: # a loop needed to obtain a single permutation
      for i in range(N-2,-2,-1): # find the largest index i where A[i] < A[i+1], if ever
        if i == -1: # no such element exists; base case
          A.reverse() # then this is the greatest permutation
          if floor(N/2.0) % 2 == 1.0: # then its odd
            parityflip = not parityflip # invert the boolean
          return even_perms # return the function
        elif A[i] < A[i+1]:
          for j in range(N-1,-1,-1): # find the largest index j such that A[i] < A[j]
            if A[i] < A[j]:
              break
          tmp = A[i] # swap A[i] and A[j]
          A[i]=A[j]
          A[j]=tmp
          parityflip = not parityflip # invert boolean
          chunk=A[i+1:] # reverse the order of A[i+1] to the end of the array
          M = len(chunk)
          chunk.reverse()
          A[i+1:]=chunk
          if floor(M/2.0) % 2 == 1.0: # then its odd
            parityflip = not parityflip
          break
      if not parityflip:
        even_perms.append(A[:])
        break
    numiter += 1
    if numiter > maxiter: # just to make sure the outer loop doesn't run out of control
      break
  return # this should never even be reached

quaternions_tesseract = [[0.5, 0.5, 0.5, 0.5],
                        [ 0.5,-0.5,-0.5,-0.5],
                        [ 0.5,-0.5, 0.5, 0.5],
                        [ 0.5, 0.5,-0.5, 0.5],
                        [ 0.5, 0.5, 0.5,-0.5],
                        [ 0.5, 0.5,-0.5,-0.5],
                        [ 0.5,-0.5,-0.5, 0.5],
                        [ 0.5,-0.5, 0.5,-0.5]]  # quaternions for the cells of a tesseract; useful for uniform SO(3) rotation

quaternions_cross_polytope = [[1, 0, 0, 0],
                              [0, 1, 0, 0],
                              [0, 0, 1, 0],
                              [0, 0, 0, 1]] # The 4D equivalent to an octahedron

quaternions_24cell = quaternions_cross_polytope + quaternions_tesseract

quaternions_simplex = normalize([[1/sqrt(10.0), 1/sqrt(6.0), 1/sqrt(3.0), 1.0],
                       [1/sqrt(10.0), 1/sqrt(6.0), 1/sqrt(3.0),-1.0],
                       [1/sqrt(10.0), 1/sqrt(6.0),-2/sqrt(3.0), 0.0],
                       [1/sqrt(10.0), -sqrt(3.0/2.0), 0.0, 0.0],
                       [-2*sqrt(2.0/5.0), 0.0, 0.0, 0.0]]) # quaternions for a pentachoron: a 4-simplex (tetrahedron in 4D)

quaternions_single = [[1.0, 0.0, 0.0, 0.0]] # quaternions for a single rotation

quaternions_test = [[1, 0, 0, 0],
                              [0.5, 0.5, 0.5, 0.5],
                              [0.5, -0.5, 0.5, 0.5]]

quaternions_tesseract1 = [] #make_SO3_grid(1)
quaternions_tesseract2 = [] #make_SO3_grid(2)

phi = (1.0 + sqrt(5.0)) / 2.0; phi_inv = 1.0 / phi
quaternions_120cell = quaternions_24cell + even_permutation([ phi/2.0, 0.5, phi_inv/2.0, 0.0]) + \
                                            even_permutation([ phi/2.0, 0.5,-phi_inv/2.0, 0.0] ) + \
                                            even_permutation([ phi/2.0,-0.5, phi_inv/2.0, 0.0] ) + \
                                            even_permutation([ phi/2.0,-0.5,-phi_inv/2.0, 0.0] )
                                            #even_permutation([-phi/2.0, 0.5, phi_inv/2.0, 0.0], reject_lower_hemisphere=True)
                                            #even_permutation([-phi/2.0, 0.5,-phi_inv/2.0, 0.0], reject_lower_hemisphere=True) # only need one side
                                            #even_permutation([-phi/2.0,-0.5, phi_inv/2.0, 0.0], reject_lower_hemisphere=True) # of the hedron
                                            #even_permutation([-phi/2.0,-0.5,-phi_inv/2.0, 0.0], reject_lower_hemisphere=True)
quaternions_120cell = quaternions_to_upper_hemisphere(quaternions_120cell)
radii = pdb.radii


def get_rot_vec(vec):
  '''takes a list of 3 numbers and makes a rotatable vector'''
  result = numpy.append(numpy.matrix([vec]).T, [[1.0]], 0)
  return result

def generate_orientation(structure,matrix):
  '''given a ligand pdb object, will rearrange the atoms in the ligand according
  to the given rotation matrix'''

  struct_center = numpy.array([pdb.center_of_mass(structure)]) # find the center of mass of a ligand
  #print struct_center
  for atom in structure.get_atoms():
    coord = get_rot_vec(atom.get_coords())
    new_coord = matrix * coord
    atom.set_coords(new_coord.T) # assign atom object to have given coordinates
  return structure

def matrix_center_of_rotation(matrix, new_center):
  '''adds a center of rotation into the rightmost column of the rotation matrix'''
  new_center = matrix(new_center)
  matrix[0:3,3] = new_center.T
  return matrix

def generate_configs(structure, locations, sites, quaternions, fullnames, traj=None,traj_matrices=None, traj_coms=None, monomer=None, random_quat=False):
  '''given a list of coordinates and a list of orientations (in the form of
    quaternions) all possible orientations at all possible locations will be
    generated'''
  if verbose: print "Generating structural configurations..."
  configs = []
  starttime = time.time()
  counter = 0
  struct_center = pdb.center_of_mass(structure)
  #for q in quaternions: # for every quaternion, generate a rotation matrix
  #  matrix = numpy.matrix(transformations.quaternion_matrix(q))
  #  matrices.append(matrix) # take the entire 4x4 matrix

  totalconfigs = len(locations) # # the total number of configs to be generated
  if verbose: print "Total number of configurations to test/generate:", totalconfigs
  #location_counter = 0
  for location in locations: # for every anchor
    q = quaternions[counter] # find the quaternion for this milestone
    if random_quat==True: # calculate a random quat
      matrix = numpy.matrix(transformations.quaternion_matrix(random_quaternion())) # generate a random quaternion inside this function
    else: # then generate a rotation matrix from the quat
      matrix = numpy.matrix(transformations.quaternion_matrix(q))
    #print "struct_center", struct_center, "location:", location
    #delta = location - struct_center
    #print "delta:",delta
    #for matrix in matrices: # for every rotation matrix

    #print "matrix:", matrix
    new_structure = deepcopy(structure)
    new_structure.moveby(-struct_center) # move the center of mass over 0,0,0
    new_structure = generate_orientation(new_structure, matrix) # rotate the structure
    #print pdb.center_of_mass(new_structure)
    new_structure.moveby(location) # move the center of mass over the specified location
    #quat_index =
    #location_index = counter
    new_structure.struct_id = fullnames[counter] #"%d_%s_%.1f_%.1f_%.1f_%d" % (counter, sites[counter], location[0],location[1],location[2], counter%len(matrices)) # give the structure a new id based on: counter, xyz-location, rotation number
    # before saving the structure, have to make sure its not clashing
    if traj: # FLAGGED FOR IMPROVEMENT
      #clashes = structures_clash_MDanalysis(new_structure, clash_structure)
      # look to see if any frame of the trajectory will allow this pose to exist
      good_frame = fit_random_snapshots(new_structure, traj, traj_range=range(5))

      if good_frame == None: # then waste no more time on this one
        print "No frame found to fit location:", location, counter
        continue
      else: # then a frame was found that will work for our purposes
        #print location, good_frame
        if verbose: print "Location", location, "fits into frame", good_frame

      new_structure.moveby(traj_coms[good_frame])
      new_structure.matrix_operation(traj_matrices[good_frame])

      #break
      if monomer == None:
        ligstring = ""
      else:
        ligstring = "%d_" % monomer
      new_structure.save('/tmp/ligand%s%0*i_frame%d.pdb' % (ligstring, int(log10(totalconfigs))+1, counter, good_frame), standard = False)
    else:
      #print "center of mass", pdb.center_of_mass(new_structure)
      #new_structure.save('/tmp/ligand%i.pdb' % counter , standard = False)
      configs.append(new_structure)
    counter += 1
    #location_counter += 1

  #print "counter: ", counter
  endtime = time.time()
  if verbose: print "Generate_configs complete. Total number %d. Elapsed time: %d" % (len(configs), endtime-starttime,)

  return configs

def structures_clash(structure1, structure2, tolerance=0.0):
  '''determines whether two structures have any clashing atoms
    the smallest structure should be the first argument to allow
    this function to run quickly as possible'''
  if verbose: print "searching for clashes between structures", structure1.struct_id, "and", structure2.struct_id
  struct1_com = pdb.center_of_mass(structure1) # find center of mass of structure1
  struct1_rad = pdb.molecular_radius(structure1) # find molecular radius
  clashing = False

  for atom2 in structure2.atoms: # for every atom in structure2
    #try:
      #atom2_radius = float(radii[atom2.resname][atom2.name])
    #except KeyError:
     # look in the pdb for the radius itself
    if atom2.radius == '0.0':
      atom2_radius = radii[atom2.element] # if not included, then retrieve a standard radius for this atom
      print "using dictionary on atom2. atom2.radius:", atom2.radius
    else:
      atom2_radius = float(atom2.radius)

    # the statement below saves some computation time by being a semi-divide and conquer method. Could be better but I'm sure it will work just fine
    if numpy.linalg.norm(numpy.array(struct1_com) - numpy.array(atom2.coords)) < struct1_rad + atom2_radius:
      # then we have to check this atom against every struct1 atom
      for atom1 in structure1.get_atoms(): # for every atom in structure1
        if atom1.radius == '0.0':
          atom1_radius = radii[atom1.element] # if not included, then retrieve a standard radius for this atom
          print "using dictionary on atom1. atom1.radius:", atom1.radius
        else:
          atom1_radius = float(atom1.radius)
        #print atom1_radius
        atom_dist = atom2_radius + atom1_radius
        if numpy.linalg.norm(numpy.array(atom1.coords) - numpy.array(atom2.coords)) < atom_dist - tolerance:
          #print numpy.linalg.norm(numpy.array(atom1.coords) - numpy.array(atom2.coords))
          # then we have a clash
          return True
  return clashing

def structures_clash_MDanalysis(structure1, MDatoms):
  '''determines whether two structures have any clashing atoms
    the smallest structure should be the first argument to allow
    this function to run quickly as possible'''
  if verbose: print "using MDAnalysis to find clash between structures", structure1.struct_id, "and", structure2.struct_id
  struct1_com = pdb.center_of_mass(structure1) # find center of mass of structure1
  struct1_rad = pdb.molecular_radius(structure1) # find molecular radius
  clashing = False
  for atom2 in MDatoms: # for every atom in structure2
    atom2_radius = radii[pdb.find_element(atom2.name)]
    if numpy.linalg.norm(numpy.array(struct1_com) - numpy.array(atom2.pos)) < struct1_rad + atom2_radius:
      # then we have to check this atom against every struct1 atom
      for atom1 in structure1.get_atoms(): # for every atom in structure1
        atom1_radius = radii[atom1.element] # radii[atom2.resname][atom2.name]
        atom_dist = atom2_radius + atom1_radius
        if numpy.linalg.norm(numpy.array(atom1.coords) - numpy.array(atom2.pos)) < atom_dist:
          return True
  return clashing

def one_eighty_matrix(axisnum):
  temp_matrix = numpy.matrix([[-1, 0,  0, 0],
                                 [ 0,-1,  0, 0],
                                 [ 0, 0, -1, 0],
                                 [ 0, 0,  0, 1]])
  temp_matrix[axisnum, axisnum] = 1
  return temp_matrix

def get_orient_matrix(struct_vec, ref_vec):
  '''used by orient_to_principal_axes function below to calculate the matrix to align two vectors'''
  rotvec = numpy.cross(struct_vec,ref_vec) # the vector around which we will rotate to align
  sine = numpy.linalg.norm(rotvec)
  cosine = numpy.dot(struct_vec,ref_vec)
  angle = atan2(sine, cosine) # the angle around the axis we will spin
  q = transformations.quaternion_about_axis(angle, rotvec)
  matrix = numpy.matrix(transformations.quaternion_matrix(q)) # convert the quaternion to a rotation matrix
  return matrix

def orient_to_axis(structure, struct_vec, ref_vec):
  '''orients the structure so that struct_vec is aligned to ref_vec'''
  #global x
  com = pdb.center_of_mass(structure)
  structure.moveby(-com)
  matrix = get_orient_matrix(struct_vec, ref_vec)
  #print "matrix:", matrix
  #x = matrix * x
  structure = generate_orientation(structure, matrix) # rotate the structure
  structure.moveby(com)
  return structure, matrix

def orient_to_principal_axes(structure):
  ''' finds the principle axes of the structure, and aligns the structure
  to the those axes'''
  evals, p_axes = pdb.principal_axes(structure) # calc principal axes
  vec_minor = p_axes[0] / numpy.linalg.norm(p_axes[0]) # the minor axis
  vec_z = numpy.array([0,0,1]) # the axis to which we will align, in this case the z
  structure, m1 = orient_to_axis(structure, vec_minor, vec_z) # orient to z
  if structure.atoms[0].coords[2] > pdb.center_of_mass(structure)[2]: # if atom0 coordinate is positive
    structure = generate_orientation(structure, one_eighty_matrix(0))
  evals, p_axes = pdb.principal_axes(structure) # recalc principal axes
  vec_inter = p_axes[1] / numpy.linalg.norm(p_axes[1]) # the intermediate axis
  vec_y = numpy.array([0,1,0]) # the axis to which we will align, in this case the y
  structure, m2 = orient_to_axis(structure, vec_inter, vec_y) # orient to y
  if structure.atoms[0].coords[1] > pdb.center_of_mass(structure)[1]: # if atom0 coordinate is positive
    structure = generate_orientation(structure, one_eighty_matrix(2))
  evals, p_axes = pdb.principal_axes(structure) # recalc principal axes
  m = m2 * m1
  #print "evals:", evals
  return m #structure, m

def random_quaternion():
  '''generates a random rotation quaternion'''
  phi = 2*pi*random()
  theta = asin(2*random() - 1)  # these assign random Euler angles of uniform orientation
  psi = pi * (2*random() - 1)
  rand_quat = [cos(phi*0.5)*cos(theta*0.5)*cos(psi*0.5) + sin(phi*0.5)*sin(theta*0.5)*sin(psi*0.5),
               sin(phi*0.5)*cos(theta*0.5)*cos(psi*0.5) - cos(phi*0.5)*sin(theta*0.5)*sin(psi*0.5),
               cos(phi*0.5)*sin(theta*0.5)*cos(psi*0.5) + sin(phi*0.5)*cos(theta*0.5)*sin(psi*0.5),
               cos(phi*0.5)*cos(theta*0.5)*sin(psi*0.5) - sin(phi*0.5)*sin(theta*0.5)*cos(psi*0.5)]
  return rand_quat

def rot_trans_test(ref_quats, test_quat):
  '''TEST FUNCTION: uses a dot product to determine which ref_quat the test_quat is closest to'''
  highest_dot = 0.0
  closest_quat = None
  counter = 0
  for ref in ref_quats:
    dot = abs(numpy.dot(numpy.array(ref),numpy.array(test_quat)))
    if dot > highest_dot:
      highest_dot = dot
      closest_quat = ref
      closest_counter = counter
    counter += 1
  #print "closest quaternion:", closest_quat, ", dot:", lowest_dot, ", counter:", closest_counter
  return closest_counter

def rot_stats(n):
  '''TEST FUNCTION'''
  quat_list = [0,0,0,0,0]
  for i in range(n):
    closest_counter = rot_trans_test(quaternions_simplex, random_quaternion())
    quat_list[closest_counter] += 1

  print "quat_list:", numpy.array(quat_list) / float(n)

def get_hedron(quat_method, quaternion_random_count=0):
  if quat_method == 'single':
    quaternions = quaternions_single
  elif quat_method == 'simplex':
    quaternions = quaternions_simplex
  elif quat_method == 'tesseract':
    quaternions = quaternions_tesseract
  elif quat_method == 'tesseract1':
    quaternions = quaternions_tesseract1
  elif quat_method == 'tesseract2':
    quaternions = quaternions_tesseract2
  elif quat_method == '24-cell':
    quaternions = quaternions_24cell
  elif quat_method == '120-cell':
    quaternions = quaternions_120cell
  elif quat_method == 'test':
    quaternions = quaternions_test
  elif quat_method == 'random':
    quaternions = []
    for i in range(quaternion_random_count): # generate a number of random quaternions
      quaternions.append(quaternions_single[0]) # because the random_quat is generated in the 'generate_configs' function when the random_quat arg is set to True
  else:
    raise Exception, "quat_method: %s not found" % quat_method
  return quaternions

def decompose_milestones(milestones):
  ''' given the compact expressions of positional and rotational milestones,
  will construct the necessary lists to make the system configs and filetree'''
  locations = []
  quaternions = []
  fullnames = []
  sites = []
  index_list = []
  positional_milestones = []
  rotational_milestones = []
  site_counter = 0
  for milestone in site:
    if milestone.shape in ["sphere","plane"]:
      positional_milestones.append(milestone)
    elif milestone.shape in ["rotational"]:
      rotational_milestones.append(milestone)
  pos_counter = 0 # counting positional milestones
  total_counter = 0 # counting all milestones
  for pos in positional_milestones:
    rot_counter = 0 # counting rotational milestones
    #if len(rotational_milestones) == 0:

    #for rot in rotational_milestones[0]:
    rot = rotational_milestones[0]
    locations.append(np.array(pos.anchor))
    quaternions.append(rot.anchor)
    sites.append(pos.siteid)
    fullname = "%d_%d_%s_%.1f_%.1f_%.1f_%d" % (total_counter, pos_counter, pos.siteid, pos.anchor[0], pos.anchor[1], pos.anchor[2], rot_counter)
    fullnames.append(fullname)
    index_list.append((pos.index, rot.index, pos.sitenum))
    rot_counter += 1
    total_counter += 1
    pos_counter += 1

  print "index_list", index_list
  return locations, quaternions, fullnames, sites, index_list

def main(settings): # NOTE: will need to include other arguments: traj, etc...
  '''parses settings and controls the rest of the script to provide a list of structures at different configurations'''
  ligand = settings['ligand']
  receptor_wet = settings['receptor']
  receptor_dry = settings['receptor_dry_pqr']
  milestones = settings['milestone_list']
  '''locations = []
  quaternions = []
  fullnames = []
  sites = []
  for m in milestones:
    locations.append(np.array(m.anchor)) # retrieve the anchor from the milestone and append it to the locations list
    quaternions.append(m.rotation)
    sites.append(m.siteid)
    fullnames.append(m.fullname)'''
  locations, quaternions, fullnames, sites, index_list = decompose_milestones(milestones)

  if settings['align_lig_to_pa']: orient_to_principal_axes(ligand)
  quat_method = settings['quaternion_method']

  if quat_method == 'random':
    random_quat = True
    #quaternions = get_hedron(quat_method, settings['quaternion_random_count'])
  else:
    random_quat = False
    #quaternions = get_hedron(quat_method)

  raw_configs = generate_configs(ligand,locations, sites, quaternions, fullnames=fullnames, random_quat=random_quat) # generate the list of structures at different configurations
  wet_configs = []
  dry_configs = []

  if verbose:
    print "Now running through all configurations to combine ligand with receptor."
    if settings['reject_clashes']:
      print "Rejecting steric clashes..."
    else:
      print "Ignoring steric clashes..."



  #included_indeces = set()
  pos_rot_combo = [] # the combinations of positions and rotations that are sterically allowed
  insert_index=0
  for i in range(len(raw_configs)): # run thru all the configurations of ligands
    lig_config = raw_configs[i]
    if verbose: print "Now creating anchor", lig_config.struct_id
    if settings['reject_clashes']:
      if structures_clash(lig_config, receptor_dry, tolerance=0.2): continue # if there's a clash
    holo_config_wet, insert_index, last_ligand_index = pdb.ligmerge(lig_config, receptor_wet, verbose=False)
    holo_config_wet.struct_id = lig_config.struct_id # set the structure description to the same as the ligand
    holo_config_wet.renumber_indeces() # to number the indeces consecutively
    wet_configs.append(holo_config_wet)
    pos_rot_combo.append(index_list[i])

  print "len(wet_configs):", len(wet_configs)
  print "len(raw_configs):", len(raw_configs)

  return wet_configs, raw_configs, pos_rot_combo, insert_index, last_ligand_index #, dry_configs

class Test_seekr(unittest.TestCase):
  # several test cases to ensure the functions in this module are working properly
  def test_normalize(self): # a small test to detect errors
    quat_list = [[1,0,0,0],[1,2,3,4]]
    expected_result = [[1,0,0,0],[0.18257418583505536, 0.3651483716701107, 0.5477225575051661, 0.7302967433402214]]
    result = normalize(quat_list)
    for i in range(len(quat_list)):
      for j in range(len(quat_list[0])):
        self.assertAlmostEqual(expected_result[i][j], result[i][j], 3)

  def test_main(self): # test this function
    print "WARNING: this module does not have comprehensive unit tests.\n\nMany of these functions require complicated inputs and outputs, that perhaps can be formulated someday."
    return
  def test_even_permutation(self):
    x = even_permutation([1,2,3,4])
    expected = [[1, 2, 3, 4], [1, 3, 4, 2], [1, 4, 2, 3], [2, 1, 4, 3], [2, 3, 1, 4], [2, 4, 3, 1], [3, 1, 2, 4], [3, 2, 4, 1], [3, 4, 1, 2], [4, 1, 3, 2], [4, 2, 1, 3], [4, 3, 2, 1]]
    self.assertEqual(x,expected)

  def test_quaternions_to_upper_hemisphere(self):
    self.assertEqual(quaternions_to_upper_hemisphere([[1,1,1,1], [2,2,2,2]]), [[1,1,1,1], [2,2,2,2]])
    self.assertEqual(quaternions_to_upper_hemisphere([[-1,-1,-1,-1], [2,2,2,2]]), [[1,1,1,1], [2,2,2,2]])
    self.assertEqual(quaternions_to_upper_hemisphere([[-1,1,-1,1], [-2,2,2,2]]), [[1,-1,1,-1], [2,-2,-2,-2]])

if __name__ == "__main__":
  print "Now running unit tests for positions_orient.py"
  unittest.main() # then run unit tests
