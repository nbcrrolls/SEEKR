'''
# hidr.py
# by Lane Votapka
# Amaro lab 2014

# HIDR (Holo Insertion by Direct Rotation): the adjunct script to be used in conjunction with SEEKR

# It takes a set of positional and rotational states. It also takes a trajectory.
# The script will run through the trajectory, and for each frame of the traj, figures out
# which state it belongs to.
#  Probably, not all states will be visited by the trajectory, and if they aren't, then
# HIDR takes the closest structure (in a different state), and exports the system configuration for
# SEEKR to place the ligand inside, run minimizations/temperature equilibrations, and then a sufficient
# ensemble equilibrations for Milestoning.


Created on Sep 18, 2014

@author: lvotapka
'''
import glob, os, sys
import xml.etree.ElementTree as ET
import MDAnalysis as mda
import numpy as np
import scipy as sp
from scipy.linalg import eigh
import unittest
import transformations # library for useful rotation functions
from random import gauss
from string import Template

k = 8.3144621454689521e-01 # in A**2 * amu / (ps**2 * K)
xsc_header = "" # this is actually pulled from the xst file

def parse_milestone_xml(milestone_filename):
  tree = ET.parse(milestone_filename)
  root = tree.getroot()
  sites = []
  N = 0 # count the number of states
  site_counter = 1
  for site in root:
    milestones = []
    for milestone in site:
      d = {}
      d['index'] = milestone.find('index').text.strip()
      d['anchor'] = np.array(map(float, milestone.find('anchor').text.strip().strip('()[]').split(',')))
      d['end'] = milestone.find('end').text.strip()
      d['bd'] = milestone.find('bd').text.strip()
      d['shape'] = milestone.find('shape').text.strip()

      d['directory'] = milestone.find('directory').text.strip()
      d['fullname'] = milestone.find('fullname').text.strip()
      if d['shape'] == 'sphere':
        d['radius'] = milestone.find('radius').text.strip()
      if d['shape'] == 'plane':
        d['normal'] = np.array(map(float, milestone.find('normal').text.strip().strip('()[]').split(' ')))
        d['distance'] = float(milestone.find('distance').text.strip())
      milestones.append(d)

      N += 1
    sites.append(milestones)
    site_counter += 1
  return sites

def confs_to_quat (coords, refcoords, com, refcom):
  n = len(coords)
  assert n == len(refcoords), "Wrong number of atoms for quaternion comparison"
  #coords = atoms.coordinates() - atoms.centerOfMass()
  #refcoords = refatoms.coordinates() - refatoms.centerOfMass()
  M = np.matrix(np.zeros((4,4)))
  for i in range(n):
    pi = coords[i] - com; qi = refcoords[i] - refcom
    Pi_trans = np.matrix([[ 0.0,   pi[0], pi[1], pi[2]],
                          [-pi[0], 0.0,  -pi[2], pi[1]],
                          [-pi[1], pi[2], 0.0,  -pi[0]],
                          [-pi[2],-pi[1], pi[0],   0.0]])
    Qi = np.matrix([[ 0,-qi[0],-qi[1],-qi[2]],
                   [ qi[0], 0,-qi[2], qi[1]],
                   [ qi[1], qi[2], 0,-qi[0]],
                   [ qi[2],-qi[1], qi[0], 0]])
    M += Pi_trans * Qi
  evecs = eigh(M)[1]
  quat = evecs[:,-1]
  return quat

def getCOM (coords, weights):
  '''given a set of coordinates and weights, will return the center of mass'''
  com = np.array([0.0,0.0,0.0])
  totalweight = 0.0
  n = len(coords)
  for i in range(n):
    coord = coords[i]
    weight = weights[i]
    combine=weight*coord
    com += combine
    totalweight += weight
  inv_totalweight = 1.0 / totalweight
  com = com*inv_totalweight
  return com

def make_random_vel(atoms, T, ts=None):
  'given a temperature and an MDAnalysis universe, will generate a velocity distribution for that temperature'
  #atoms = u.selectAtoms('all')
  n = len(atoms)
  velcoords = np.zeros((n,3))
  for i in range(n):
    mass = atoms[i].mass
    velcoords[i,0] = gauss(0.0, np.sqrt(k*T / mass))
    velcoords[i,1] = gauss(0.0, np.sqrt(k*T / mass))
    velcoords[i,2] = gauss(0.0, np.sqrt(k*T / mass))
  atoms.positions = velcoords
  return

def extract_xst_frames(xst_filename_list, filename_suffix_list, directory_list, rootdir, xsc_prefix, start=0, end=-1, stride=1):
  'extracts all xsc files from an series of xst trajs that should correspond to dcd trajs'

  total_index = 0
  i = -3
  for xst_filename in xst_filename_list:
    xst_file = open(xst_filename,'r')
    file_index = -3
    for line in xst_file.xreadlines():
      if file_index > 0 and i >= start and (i < end or end == -1) and (i-start)%stride == 0:
        dirname = directory_list[i]
        xsc_filename = os.path.join(rootdir, dirname, "md/reverse", "%s%s.xsc" % (xsc_prefix, filename_suffix_list[i]))
        xsc_file = open(xsc_filename, 'w')
        xsc_file.write(xsc_header)
        xsc_file.write(line)
        xsc_file.close()
        total_index += 1
      elif i == -2:
        xsc_header = line
      elif i == -1:
        xsc_header += line
      i += 1
      file_index += 1

class Test_hidr(unittest.TestCase):
  # several test cases to ensure the functions in this module are working properly
  def test_main(self): # test this function
    print "WARNING: this module does not have comprehensive unit tests"

  def test_getCOM(self):
    self.assertTrue(np.allclose(getCOM(np.array([[0.5, 0.0, 2.5],[3.5, 0.0, 2.5],[2.0, 0.0, 0.5]]), [3.0, 3.0, 4.0]), np.array([2.0, 0.0, 1.7])))
    self.assertTrue(np.allclose(getCOM(np.array([[0.5, 2.5, 0.0],[3.5, 2.5, 0.0],[2.0, 0.5, 0.0]]), [3.0, 3.0, 4.0]), np.array([2.0, 1.7, 0.0])))

  def test_confs_to_quat(self):
    atoms=np.array([[1, 1, -1], [2, -3, 4], [5, -6, 7], [-5, 7, 0], [-5, 6, 5]])
    refatoms=np.array([[1, 1, -1],[2, -3, 4],[5, -6, 7],[-5, 7, 0],[-5, 6, 5]])
    weights=[2.0, 2.7, 1.0, 1.2, 0.5]
    com=getCOM(atoms, weights)
    refcom=getCOM(refatoms, weights)
    quat=confs_to_quat(atoms, refatoms, com, refcom)
    self.assertTrue(np.allclose(quat , np.array([1.0, 0.0, 0.0, 0.0])), msg="quaternions are not similar") # no rotation
    atoms=np.array([[-2.013885498046875, -2.112818956375122, 0.7542771100997925],[4.031986236572266, -0.8651127815246582, 2.7267541885375977],[9.1787748336792, -0.6724526882171631, 2.0386807918548584],[-8.197371482849121, 2.896355390548706, 3.8643691539764404],[-4.601040840148926, 6.116539001464844, 5.506570816040039]])
    com=getCOM(atoms, weights)
    quat = confs_to_quat(atoms, refatoms, com, refcom)
    self.assertTrue(np.allclose(quat, np.array([0.7662610281769211, 0.38313051408846055, -0.47891314261057566, -0.19156525704423027])), msg="quaternions are not similar") # random rotation with same COM
    atoms=np.array([[-6.513885498046875, -8.812818956375121, 3.0542771100997923],[-0.4680137634277344, -7.565112781524658, 5.0267541885375975],[4.678774833679199, -7.372452688217163, 4.338680791854858],[-12.697371482849121, -3.803644609451294, 6.16436915397644],[-9.101040840148926, -0.5834609985351564, 7.806570816040039]])
    com=getCOM(atoms, weights)
    self.assertTrue(np.allclose(confs_to_quat(atoms, refatoms, com, refcom), np.array([0.7662610281769211, 0.38313051408846055, -0.47891314261057566, -0.19156525704423027])), msg="quaternions are not similar") ;# random rotation with different COM
    return

#if __name__ == '__main__':
# Define all the necessary files to load the long trajectory
# the directory where the SEEKR project resides
root_dir="/scratch/lvotapka/projects/seekr/urea_membrane_rot"
delete_old_namd_files = True

# the trajectory files to load
traj_pdb="/extra/moana/amarolab1/Permeability_Projects/simulations/Systems/urea/urea_dmpc.pdb"
traj_psf="/extra/moana/amarolab1/Permeability_Projects/simulations/Systems/urea/urea_dmpc.psf"
#set traj_list { \
#"/extra/moana/amarolab1/Permeability_Projects/simulations/Production/urea/window_-35/umb-1.dcd" \
#}
#traj_list=glob.glob('/extra/moana/amarolab1/Permeability_Projects/simulations/Production/urea/*/umb-?.dcd') # find all trajectories
traj_list=glob.glob('/scratch/lvotapka/projects/seekr/urea_membrane/anchor_*/md/ens_equil/ens_equil_*.dcd') # NOTE: change all these to the fixed dcds in my scratch directory
veltraj_list=glob.glob('/scratch/lvotapka/projects/seekr/urea_membrane/anchor_*/md/ens_equil/ens_equil_*.veldcd') # NOTE: change all these to the fixed dcds in my scratch directory
xst_filename_list = glob.glob('/scratch/lvotapka/projects/seekr/urea_membrane/anchor_*/md/ens_equil/ens_equil_*.xst')

ligstr = 'resname URE'
recstr = 'name P'
# load the reference structures
ignore_rec_rot = True # ignore the receptor's rotational states (say, for a membrane)



# define the states to compare the traj to
state_dict = {}
# READ MILESTONES.XML FILE
rotational_milestones = []
positional_milestones = []
milestone_xml_filename = os.path.join(root_dir,'milestones.xml')
sites = parse_milestone_xml(milestone_xml_filename)
for milestone_list in sites: # each site has its own set of milestones
  for milestone in milestone_list: # sort each milestone as to whether it's rotational/positional
    if milestone['shape'] == "rotational":
      rotational_milestones.append(milestone)
    else:
      positional_milestones.append(milestone)

for pos in positional_milestones: # iterate thru all possible combinations of milestone positions and rotations
  com = pos['anchor']
  for rot in rotational_milestones:
    quat = rot['anchor']
    key = "_".join((str(com),str(quat))) # fill the state dictionary's keys
    state_dict[key] = []
# make a list of .xst files from traj_list
for trajfile in traj_list:
  basename = os.path.splitext(trajfile)[0] # the file without the extension

#print "positional_milestones:", positional_milestones

# open the trajectory
md_universe = mda.Universe(traj_psf, traj_list)
ref_universe = mda.Universe(traj_psf, traj_pdb)
ligand = md_universe.selectAtoms(ligstr)
ligref = ref_universe.selectAtoms(ligstr) # the reference ligand orientation
receptor = md_universe.selectAtoms(recstr)
recref = ref_universe.selectAtoms(recstr)
lig_ref_quat = confs_to_quat(ligand.atoms.coordinates(), ligref.atoms.coordinates(), ligand.atoms.centerOfMass(), ligref.atoms.centerOfMass())
rec_ref_quat = confs_to_quat(receptor.atoms.coordinates(), recref.atoms.coordinates(), receptor.atoms.centerOfMass(), recref.atoms.centerOfMass())
pos_all = md_universe.selectAtoms("all") # select all the atoms to write to the position file
vel_all = md_universe.selectAtoms("all") # select all the atoms to write to the velocity file (random velocities)


#md_universe
counter = 0
indexpair_list = [] # a list of the rotation-counter indeces
directory_list = [] # list of the directories of each index
directory_set = set() # a set of the directories that we are going thru
for ts in md_universe.trajectory:
  # find the COM & orientation of the ligand and receptor
  ligCOM = ligand.centerOfMass()
  recCOM = receptor.centerOfMass()
  lig_quat = confs_to_quat(ligand.atoms.coordinates(), ligref.atoms.coordinates(), ligand.atoms.centerOfMass(), ligref.atoms.centerOfMass())
  rec_quat = confs_to_quat(receptor.atoms.coordinates(), recref.atoms.coordinates(), receptor.atoms.centerOfMass(), recref.atoms.centerOfMass())

  # find which state it belongs to
  lowest_dist = sys.maxint
  pos_state = None

  for pos in positional_milestones:
    milestoneCOM = recCOM + pos['normal']*pos['distance']
    rel_coord = ligCOM - milestoneCOM
    if pos['shape'] == 'plane':
      dist = np.dot(rel_coord, pos['normal'])
    if pos['shape'] == "sphere":
      radius = float(pos['radius'])
      dist = np.dot(rel_coord, pos['normal']) - radius
    # see if this is the closest milestone so far
    if abs(dist) < lowest_dist: # then we've found a closer milestone
      lowest_dist = abs(dist)
      pos_state = pos

  #print "pos_state:", pos_state['anchor']

  highest_similarity = 0.0
  highest_index = -1
  rot_state = None
  #rot_index = 0
  if ignore_rec_rot: rec_quat = [1,0,0,0]
  for rot in rotational_milestones:
    # Find the quaternion representing this timestep's rotation
    timestep_quat = transformations.quaternion_multiply(lig_quat, transformations.quaternion_inverse(rec_quat))
    # Find the quaternion representing this rotational milestone, compare to the timestep quat
    rot_index = int(rot['index'])
    similarity = np.dot(timestep_quat, rot['anchor'])
    if abs(similarity) > highest_similarity: # then we've found a closer milestone
      highest_similarity = abs(similarity)
      rot_state = rot
      highest_index = rot_index
    #rot_index += 1
  #print "rot_state:", rot_state['anchor']
  key = "_".join((str(pos_state['anchor']),str(rot_state['anchor']))) # fill the state dictionary's keys
  state_dict[key].append(counter)

  # write the pos file
  directory = pos_state['directory']
  pos_all.write(os.path.join(root_dir, directory, "md/reverse", "pos%d_%d.pdb" % (highest_index, counter)))
  make_random_vel(vel_all,counter)
  vel_all.write(os.path.join(root_dir, directory, "md/reverse", "vel%d_%d.pdb" % (highest_index, counter)))
  indexpair_list.append("%d_%d" % (highest_index, counter))
  directory_list.append(directory)
  directory_set.add(directory)
  counter += 1

# create a file that lists all the anchors that need more processing
for key in state_dict.keys():
  if len(state_dict[key]) == 0:
    print "Alert: no frames found for state:", key

# make the .namd files
counter = 0
old_dir_filename = ""

for dir_filename in directory_set:
  if delete_old_namd_files:
    namd_files_to_delete = glob.glob(os.path.join(root_dir,dir_filename,"md/reverse","reverse*.namd")) # Maybe instead of this, we could leave the .namd files that don't have a state yet
    #print "Deleting these files:", namd_files_to_delete
    for namd_file in namd_files_to_delete:
      os.remove(namd_file)

# load the template file
for index_pair in indexpair_list: # NOTE: what am I going to do about all these other ones here? Maybe have a hidr=True option that won't generate the namd files in the reversals
  dir_filename = directory_list[counter]
  if dir_filename != old_dir_filename:
    old_dir_filename = dir_filename
    namd_template_filename = (os.path.join(root_dir, dir_filename, 'md/reverse', "reverse_template"))
    namd_template_file = open(namd_template_filename,'r')
    namd_template_string=''.join(namd_template_file.readlines())
    namd_template_file.close()
    namd_template = Template(namd_template_string)

  # fill all the values in the template file
  index_pair_list = index_pair.split('_')
  namd = namd_template.safe_substitute(ROTINDEX=index_pair_list[0], TRAJID=index_pair_list[1])
  # write the namd file
  namd_outname = os.path.join(root_dir, dir_filename, 'md/reverse', "reverse%s.namd" % index_pair)
  namd_out = open(namd_outname, 'w')
  namd_out.write(namd)
  namd_out.close()

  counter += 1

# make the xsc files
'''
print "xst_filename_list:", xst_filename_list
print "len(xst_filename_list):", len(xst_filename_list)
print "indexpair_list:", indexpair_list
print "len(indexpair_list):", len(indexpair_list)
print "directory_list:", directory_list
print "len(directory_list):", len(directory_list)
'''
extract_xst_frames(xst_filename_list, indexpair_list, directory_list, root_dir, "reverse", start=0, end=len(directory_list), stride=1)

# now all existing states are correctly created! need to figure out which states are missing and then code in the targetedMD stuff to get those going


