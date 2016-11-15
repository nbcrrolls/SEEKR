#!/usr/bin/python

'''
seekr.py is the control program for the SEEKR package

SEEKR (Simulation Enabled Estimation of Kinetic Rates) is a package of scripts that prepare and analyze simulations for the
estimation of kinetic rates (and thermodynamics) of chemical processes.

'''

import pdb2 as pdb
import milestones
import positions_orient
import filetree
import md
import bd
import numpy as np
import unittest
import colvars
import os, shutil, sys, re
#import pickle
import cPickle as pickle
import argparse

########################################################################################
# define all SEEKR parameters passed to other programs
########################################################################################


def boolean(arg):
  if str(arg).lower() in ("false","", "0", "0.0", "[]", "()", "{}"):
    return False
  else:
    return True

#if len(sys.argv) <= 1:
#  unittest.main()
#  exit()


md_time_factor = 2.0 # 2 fs per timestep
bd_time_factor = 1000.0 # 1000 fs per ps

inp = { # contains default parameters in case they aren't included in the input file
  'test_mode':False, # reduces the calculation time by restricting input size. Use only for debugging SEEKR
  'LA_src':"../../la.tcl",
  'tcl_script_control':True, # whether the milestoning is controlled by the TCL script versus something else
  #'planescale':'1.0',
  'rec_rot':True,
  'recrot':"True",
  'lig_rot':False,
  'ligrot:':'False',
  'md':True,
  'bd':True,
  'k_off':False,
  'receptor_type':'globular',
  'empty_rootdir':False,
  'extract_stride':'1',
  'rec_xsc_filename':'',
  'script_interval':'5',
  'abort_on_crossing':False,
  #'ligpa1_list':"", # MARKED FOR REMOVAL
  #'ligpa3_list':"", # MARKED FOR REMOVAL
  'align_lig_to_pa':'False',
  'lig_psf_filename':'',
  'hedron':'single',
  'watermodel':'',
  'cell_shape': 'box',
  'quaternion_random_count':'1',
  'reject_clashes':"False",
  'quat_lig_indeces':"",
  'quat_rec_indeces':"",
  'quat_force_constant':"0.0",
  'tclforcesscript':'../../../milestoning_new.tcl',
  'quatforcesscript':'../../../quat_forces.tcl',
  'recrange':"",
  'rec_com_indeces':"",
  'rec_psf_filename':"",
  'quatforcesscript':'',
  #'recpa1_list':"", # MARKED FOR REMOVAL
  #'recpa3_list':"", # MARKED FOR REMOVAL
  'include_only_adjacent_milestones':True,
  'leap_preload_commands':'',
  'leap_postload_commands':'',
  'leap_program':'tleap',
  'sample_leap_file':'',
  'charmm_parameters':"[]",
  'recpsf':'',
  'ligpsf':'',
  'ignore_psf':False,
  'min_constrained':[],
  'min_restrained':[],
  'min_restrained_force':'10.0',
  'temp_equil_constrained':[],
  'temp_equil_restrained':[],
  'temp_equil_restrained_force':'10.0',
  'ens_equil_constrained':[],
  'ens_equil_restrained':[],
  'ens_equil_restrained_force':'10.0',
  #'ens_equil_colvar_type':'sphere',
  'ens_equil_colvar_sel': "[ 'ligand', 'receptor' ]",
  'ens_equil_colvars':True,
  'empty_pqrxml_path': "./emtpy.pqrxml",
  'fhpd_numtraj':'1000',
  'extract_xst':"True",
  'empty_pqrxml_path':"./empty.pqrxml",
  'fwd_rev_constrained':[],
  'fwd_rev_max_num_steps':'0',
  'fwd_rev_launches_per_config':'1',
  'fwd_rev_frame_chunk_size':'1000',
  'browndye_bin_dir':'',
  'ion1rad':None,
  'ion2rad':None,
  'ion1conc':None,
  'ion2conc':None,
  'ions':[],
  'lpbe_npbe':'lpbe',
  'bd_rec_pqr_filename':'',
  'apbs_executable':'apbs',
  'inputgen_executable':'inputgen.py',
  'inputgen_fadd':'100',
  'inputgen_gmemceil':'64000',
  'inputgen_resolution':'0.5',
  'inputgen_cfac':'4.0',

}

def pickle_or_load(filename, picklename, struc_name="pickle",pqr=False):
  'for large files, instead of parsing, they can be saved and loaded much more quickly as a pickle. '
  if os.path.exists(picklename) and os.path.getmtime(picklename) > os.path.getmtime(filename): # if the pickle has been most recently modified
    # load the pickle
    print "reading pickle:", picklename
    our_file=open(picklename, 'rb')
    our_obj=pickle.load(our_file)
    our_file.close()
  else:
    # then load the file itself and save the pickle
    our_obj=parser.get_structure(struc_name, filename, pqr=pqr, conventional=False) # load the file
    print "writing pickle:", picklename
    our_file=open(picklename, 'wb')
    pickle.dump(our_obj, our_file, protocol=-1)
    our_file.close()
  return our_obj

# Parse the SEEKR input file
def parse_input_file(inp_filename):
  new_inp = {}
  listopen = False # this gives us the ability to read parameters across several lines
  ourlist = []
  for line in open(inputscript,'r'):
    splitline = line.strip().split() + ['#']

    if not listopen: splitlist = []
    for i in range(len(splitline)):
      if splitline[i][0] == '[':
        splitline[i] = splitline[i][1:].strip()
        listopen = True
        if splitline[i] and splitline[i][-1] == ']':
          listopen = False
          break
      if splitline[i] and splitline[i][-1] == ']':
        listopen = False
        splitlist.append(splitline[i][:-1])
        break

      if splitline[i].startswith('#'): # the remainder of this line is a comment
        break
      splitlist.append(splitline[i].strip())

    if not listopen:
      if len(splitlist) < 2: continue
      word1 = splitlist[0].lower()
      word2 = " ".join(splitlist[1:])
      if word1[0] == '#': # this line is a comment
        continue
      new_inp[word1] = word2 # populate the input dictionary with this value
  return new_inp

parser = argparse.ArgumentParser(description="The program that is called to prepare MD/BD for a ligand-receptor system. The first argument defines the input file. If the second argument is '-r' or '--remove', it will rewrite an entire SEEKR filetree: similar to the 'empty_rootdir' option for the input file defined below. If seekr.py is run without arguments, the unit tests will be performed.")
parser.add_argument('input', metavar='INPUT', type=str, nargs='?', help="The SEEKR input file containing all the options, inputs, and parameters for the SEEKR calculation.")
# optional arguments
remove_group = parser.add_mutually_exclusive_group()
remove_group.add_argument('-r', '--remove', dest='remove', default=False, help="Remove the currently existing SEEKR file tree to be rewritten by the new files.", action="store_true")
  #parser.add_argument('-r', '--remove', metavar="REMOVE", dest='remove', default='', type=str, help='The time to run the simulations.')
args = parser.parse_args()
args = vars(args) # convert to a dictionary  

# combine the default input with the user-defined input
#inputscript = sys.argv[1] # the input file for SEEKR
if args['input']:
  inputscript = args['input']
  new_inp = parse_input_file(inputscript)
  inp.update(new_inp)
else:
  unittest.main()
  exit()

# Start defining necessary parameters from the input file
test_mode=boolean(inp['test_mode']) # reduces the calculation time by restricting input size
master_temperature = float(inp['master_temperature'])
ens_equil_len = int(inp['ens_equil_len']) # number of timesteps
number_of_ens_equil_frames = int(inp['number_of_ens_equil_frames']) # number of frames to write after the ens_equil_simulations
number_of_ens_equil_frames_skipped = int(inp['number_of_ens_equil_frames_skipped'])
extract_stride = int(inp['extract_stride'])
#number_of_ens_equil_frames_skipped = number_of_ens_equil_frames = number_of_reversals # number of frames to skip
dcd_freq = ens_equil_len / number_of_ens_equil_frames # the frequency to write dcd files in ens_equil_simulation

MILESTONE_FILENAME = os.path.join(inp['rootdir'], "milestones.xml")

inp['ligrot'] = inp['lig_rot']
inp['recrot'] = inp['rec_rot']

sys_params={ # variables pertinent to the receptor/ligand
  'project_name':inp['project_name'],
  'rootdir':inp['rootdir'],
  'lig_pdb_filename':inp['lig_pdb_filename'],
  'lig_pqr_filename':inp['lig_pqr_filename'],
  'rec_pdb_filename':inp['rec_pdb_filename'],
  'rec_psf_filename':inp['rec_psf_filename'],
  'rec_dry_pdb_filename':inp['rec_dry_pdb_filename'],
  'rec_dry_pqr_filename':inp['rec_dry_pqr_filename'],

  'rec_rot':boolean(inp['recrot']), # True: the milestones will rotate with the receptor structure. False: the milestones will remain where they are in space
  'lig_rot':boolean(inp['ligrot']), # NOTE: boolean form of this variable. True: ligand rotation transitions will be counted. False: Only translational transitions of the ligand will be noted NOTE: not yet implemented
  'md':boolean(inp['md']),
  'bd':boolean(inp['bd']),
  'bd_rec_pqr_filename':inp['bd_rec_pqr_filename'],
  'empty_rootdir':boolean(inp['empty_rootdir']), # if set to True, will empty the contents of the rootdir when the new file tree is made

}

#define program path variables

#program_paths={# variables that define paths to programs used by SEEKR
#  'namd_special':inp['namd_special'],
#  'charm_special':inp['charm_special'],
#  'mpiexec':inp['mpiexec'],

#}

#if 'remove' in sys.argv[2:]: 
if args['remove'] == True: # then remove the entire directory as if empty_rootdir was True
  sys_params['empty_rootdir'] = True

#if "milestone_xml" in sys.argv[2:]: # then skip a bunch of steps straight to making the milestoning file

print "Creating file tree at location: %s" % inp['rootdir']
if sys_params['empty_rootdir'] and os.path.isdir(inp['rootdir']): # then we're in test mode, delete anything that is in the rootdir
  print "'empty_rootdir' set to True. Deleting all subdirectories/files in rootdir:", inp['rootdir']
  shutil.rmtree(inp['rootdir']) # remove all the files within the grid subdirectory
if not os.path.isdir(inp['rootdir']): # if the rootdir doesn't exist yet
  os.mkdir(inp['rootdir']) # then create it

# find the recrange and ligrange values


if md:
  tcl={ # all variables that will go into the tcl script
    'script_interval':inp['script_interval'],
    'abort_on_crossing':inp['abort_on_crossing'],
    'ligrange':inp['ligrange'], # what actually defines the ligand
    'lig_com_indeces':inp['lig_com_indeces'], # the residues we monitor for center of mass in the simulation
    'quat_lig_indeces':inp['quat_lig_indeces'],
    'recrot':inp['recrot'],
    'ligrot':inp['ligrot'],
    #'ligpa1_list':inp['ligpa1_list'], # MARKED FOR REMOVAL
    #'ligpa3_list':inp['ligpa3_list'], # MARKED FOR REMOVAL
    'recrange':inp['recrange'],
    'rec_com_indeces':inp['rec_com_indeces'],
    'quat_rec_indeces':inp['quat_rec_indeces'],
    'quat_force_constant':inp['quat_force_constant'],
    'quatforcesscript':inp['quatforcesscript'],
    #'recpa1_list':inp['recpa1_list'], # MARKED FOR REMOVAL
    #'recpa3_list':inp['recpa3_list'], # MARKED FOR REMOVAL
  }
else:
  tcl={}

######### Milestones ##########################################################


# find all sites
sites = []
for key in sorted(inp.keys()):
  if re.match("site[0-9]+", key): # then this is a site
    site_dict = {'key':key}
    for param in inp[key].split(','): # then pull every string in the list and make a
      split_param = param.strip().split()
      if len(split_param) < 2: continue
      site_key = split_param[0]
      site_value = ' '.join(split_param[1:])
      site_dict[site_key] = site_value
    sites.append(site_dict)


# classify and log sites
milestone_settings={'sites':[]}
absolute_mode = False
for site_dict in sites:
  absolute = "False"
  startvx = 0.0; startvy = 0.0; startvz = 0.0;
  if 'startvx' in site_dict: startvx = float(site_dict['startvx'])
  if 'startvy' in site_dict: startvy = float(site_dict['startvy'])
  if 'startvz' in site_dict: startvz = float(site_dict['startvz'])
  if 'radius_list' in site_dict: 
    radius_string= site_dict['radius_list']
    raw_radius_list= radius_string.split()
    print 'radius list 1', raw_radius_list 
    radius_list=map(float,raw_radius_list)
  else: radius_list=None
  print 'radius list' , radius_list   
  if 'hedron' in site_dict: # this is so the user can easily leave out this argument, SEEKR will take it from elsewhere in the input
    milestone_hedron = site_dict['hedron']
  else:
    milestone_hedron = inp['hedron'] # straight from the input file

  if site_dict['anchor_function'] == 'concentric_spheres_atom':
    if 'absolute' in site_dict.keys() and site_dict['absolute'].lower() == 'true':
      absolute="True"
      absolute_mode = True
    milestone_settings['sites'].append({ # parameters pertaining to the milestones
      'anchor_function':site_dict['anchor_function'],
      'dimensions':{'r':float(site_dict['r']),
                  'r_low':float(site_dict['r_low']),
                  'x':float(site_dict['x']),
                  'y':float(site_dict['y']), # sphere center, in this case 'concentric_spheres_atom', only defines the starting locations of the trajectories
                  'z':float(site_dict['z']),
                  'atomid':"'%s'" % site_dict['atomid'],
                  'vx':float(site_dict['vx']),
                  'vy':float(site_dict['vy']),
                  'vz':float(site_dict['vz']),
                  'startvx':startvx,
                  'startvy':startvy,
                  'startvz':startvz,
                  'increment':float(site_dict['increment']), # increment between anchors
                  'siteid':"'%s'" % site_dict['key'], # have to put quotes otherwise milestones.py won't recognize it as a string
                  'absolute':absolute,
                  'k_off':boolean(inp['k_off']),
                  'radius_list':radius_list,
                  },

    })

  if site_dict['anchor_function'] == 'concentric_spheres_atom_with_rotations':
    if 'absolute' in site_dict.keys() and site_dict['absolute'].lower() == 'true':
      absolute="True"
      absolute_mode = True
    milestone_settings['sites'].append({ # parameters pertaining to the milestones
      'anchor_function':site_dict['anchor_function'],
      'dimensions':{'r':float(site_dict['r']),
                  'r_low':float(site_dict['r_low']),
                  'x':float(site_dict['x']),
                  'y':float(site_dict['y']), # sphere center, in this case 'concentric_spheres_atom', only defines the starting locations of the trajectories
                  'z':float(site_dict['z']),
                  'atomid':"'%s'" % site_dict['atomid'],
                  'vx':float(site_dict['vx']),
                  'vy':float(site_dict['vy']),
                  'vz':float(site_dict['vz']),
                  'startvx':startvx,
                  'startvy':startvy,
                  'startvz':startvz,
                  'increment':float(site_dict['increment']), # increment between anchors
                  'siteid':"'%s'" % site_dict['key'], # have to put quotes otherwise milestones.py won't recognize it as a string
                  'absolute':absolute,
                  'k_off':boolean(inp['k_off']),
                  'hedron':"'%s'" % milestone_hedron,
                  'radius_list':radius_list,
                  },

    })

  if site_dict['anchor_function'] == 'concentric_spheres':
    if 'absolute' in site_dict.keys() and site_dict['absolute'].lower() == 'true':
      absolute="True"
      absolute_mode = True
    milestone_settings['sites'].append({ # parameters pertaining to the milestones
      'anchor_function':site_dict['anchor_function'],
      'dimensions':{'r':float(site_dict['r']),
                  'x':float(site_dict['x']),
                  'y':float(site_dict['y']), # sphere center, in this case 'concentric_spheres_atom', only defines the starting locations of the trajectories
                  'z':float(site_dict['z']),
                  'r_low':float(site_dict['r_low']),
                  'vx':float(site_dict['vx']),
                  'vy':float(site_dict['vy']),
                  'vz':float(site_dict['vz']),
                  'startvx':startvx,
                  'startvy':startvy,
                  'startvz':startvz,
                  'increment':float(site_dict['increment']), # increment between anchors
                  'siteid':"'%s'" % site_dict['key'], # have to put quotes otherwise milestones.py won't recognize it as a string
                  'absolute':absolute,
                  'k_off':boolean(inp['k_off']),
                  'radius_list':radius_list,
                  },

    })
  if site_dict['anchor_function'] == 'planes':
    if 'absolute' in site_dict.keys() and site_dict['absolute'].lower() == 'true':
      absolute="True"
      absolute_mode = True
    milestone_settings['sites'].append({ # parameters pertaining to the milestones
      'anchor_function':site_dict['anchor_function'],
      'dimensions':{
                  'origx':float(site_dict['origx']),
                  'origy':float(site_dict['origy']), # sphere center, in this case 'concentric_spheres_atom', only defines the starting locations of the trajectories
                  'origz':float(site_dict['origz']),
                  'normx':float(site_dict['normx']),
                  'normy':float(site_dict['normy']),
                  'normz':float(site_dict['normz']),
                  'lowest':float(site_dict['lowest']),
                  'highest':float(site_dict['highest']),
                  'increment':float(site_dict['increment']), # increment between anchors
                  'siteid':"'%s'" % site_dict['key'], # have to put quotes otherwise milestones.py won't recognize it as a string
                  'absolute':absolute,
                  'k_off':boolean(inp['k_off']),
                  'radius_list':radius_list,
                  },

    })
  if site_dict['anchor_function'] == 'planes_with_rotations':
    if 'absolute' in site_dict.keys() and site_dict['absolute'].lower() == 'true':
      absolute="True"
      absolute_mode = True
    milestone_settings['sites'].append({ # parameters pertaining to the milestones
      'anchor_function':site_dict['anchor_function'],
      'dimensions':{
                  'origx':float(site_dict['origx']),
                  'origy':float(site_dict['origy']), # sphere center, in this case 'concentric_spheres_atom', only defines the starting locations of the trajectories
                  'origz':float(site_dict['origz']),
                  'normx':float(site_dict['normx']),
                  'normy':float(site_dict['normy']),
                  'normz':float(site_dict['normz']),
                  'lowest':float(site_dict['lowest']),
                  'highest':float(site_dict['highest']),
                  'increment':float(site_dict['increment']), # increment between anchors
                  'siteid':"'%s'" % site_dict['key'], # have to put quotes otherwise milestones.py won't recognize it as a string
                  'absolute':absolute,
                  'k_off':boolean(inp['k_off']),
                  'hedron':"'%s'" % milestone_hedron,
                  'radius_list':radius_list,
                  },

    })
  elif site_dict['anchor_function'] == 'ellipsoid':
    pass # NOTE: needs to be filled out later
milestone_settings['milestone_filename'] = MILESTONE_FILENAME
milestone_settings.update({'include_only_adjacent_milestones':inp['include_only_adjacent_milestones'], # only include adjacent milestone in every NAMD TCL script. This will be set to True if we will stop the simulation as soon as another milestone is reached. Set to False if simulations will not be stopped once a new milestone is reached, then every single milestone/anchor pair will need to be included in the NAMD TCL script
  'test_mode':test_mode,})



milestone_settings['master_temperature'] = master_temperature

raw_milestone_list = milestones.main(milestone_settings)
print "len(raw_milestone_list): ",len(raw_milestone_list)
#milestone_settings['milestone_list'] = milestone_list
######### Structures ############################################################
parser = pdb.Big_PDBParser()
print "now loading structures"

# Read and/or create pickle files for the structures to save I/O time
ligand_pkl_filename = os.path.join(inp['rootdir'], "ligand.pkl")
receptor_pkl_filename = os.path.join(inp['rootdir'], "receptor.pkl")
receptor_pkl_dry_filename = os.path.join(inp['rootdir'], "receptor_dry.pkl")
receptor_pkl_dry_pqr_filename = os.path.join(inp['rootdir'], "receptor_dry_pqr.pkl")

ligand=pickle_or_load(sys_params['lig_pqr_filename'], ligand_pkl_filename, struc_name="ligand", pqr=True)
receptor=pickle_or_load(sys_params['rec_pdb_filename'], receptor_pkl_filename, struc_name="receptor", pqr=False)
receptor_dry=pickle_or_load(sys_params['rec_dry_pdb_filename'], receptor_pkl_dry_filename, struc_name="receptor_dry", pqr=False)
receptor_dry_pqr=pickle_or_load(sys_params['rec_dry_pqr_filename'], receptor_pkl_dry_pqr_filename, struc_name="receptor_dry_pqr", pqr=True)

struct={ # all parameters pertaining to structure
  'ligand':ligand,
  'receptor':receptor,
  'receptor_dry':receptor_dry, # or create a function that will remove all waters, complicated by ions
  'receptor_dry_pqr':receptor_dry_pqr,
  'rec_com':pdb.center_of_mass(receptor_dry), # have to take into account the center of mass of the receptor itself
  'lig_com':pdb.center_of_mass(ligand), # have to take into account the center of mass of the ligand itself

}
print "rec_com:", struct['rec_com'], ", lig_com:", struct['lig_com']
print "ligand index", ligand.atoms[-1].index
######## Positions/Orientations ###################################################
pos_settings={ # position/orientation settings
  'quaternion_method':inp['hedron'], # Options: single, random, simplex, tesseract, ... more to come someday
  'quaternion_random_count':int(inp['quaternion_random_count']), # if quaternion_method is set to 'random', must include number of random rotations to generate per position
  'align_lig_to_pa':boolean(inp['align_lig_to_pa']), # whether the ligand structure should be aligned to principal axes (True), or left as is (False)
  'reject_clashes':boolean(inp['reject_clashes']), # True: run a check to see whether ligand/receptor structures are clashing, False: disregard clashes

}
if test_mode:
  raw_milestone_list=raw_milestone_list[-1:] # reduces input size
  pos_settings['quaternion_method'] = inp['hedron'] # only one rotation, just to make things faster

pos_settings_all = dict(pos_settings.items() + struct.items() + [('milestone_list',raw_milestone_list)])
wet_configs, lig_configs, pos_rot_combo, insert_index, last_ligand_index = positions_orient.main(pos_settings_all) # wet_configs: ligand+receptor structures, lig_configs: ligand configs only

if tcl['ligrange'].lower() == 'auto':
  # then we have to find the actual indeces
  tcl['ligrange'] = "%d to %d" % (insert_index+1, insert_index+1+last_ligand_index)
  #print "REPLACING ligrange with:", tcl['ligrange']


if tcl['recrange'].lower() == 'auto':
  # then we have to find the actual indeces
  tcl['recrange'] = "%d to %d" % (1, insert_index)
  #print "REPLACING recgrange with:", tcl['recrange']


if tcl['lig_com_indeces'].lower() == 'auto':
  #tcl['lig_com_indeces'] = tcl['ligrange']
  tcl['lig_com_indeces'] = "%d to %d" % (insert_index+1, insert_index+1+last_ligand_index)


if tcl['rec_com_indeces'].lower() == 'auto_all':
  #tcl['rec_com_indeces'] = tcl['recrange']
  tcl['rec_com_indeces'] = "%d to %d" % (1, insert_index)

if tcl['rec_com_indeces'].lower() == 'auto_ca':
  ca_indeces = [] # indeces for the alpha carbons
  for atom in wet_configs[0].atoms: # by default use the first wet structure
    if atom.name == "CA": # then it's an alpha carbon
      ca_indeces.append(atom.index)

  tcl['rec_com_indeces'] = ' '.join(map(str,ca_indeces))
  #print "REPLACING rec_com_indeces with:", tcl['rec_com_indeces']

# need to filter out the members of milestone_list so we are only including the members that have been included in the configs
print "position/rotation index combinations:", pos_rot_combo
milestone_pos_rot_list = [] # contains pairs of milestone objects for position and rotation
#milestone_pos_rot_list_indeces = [] # same as above, except contains the milestone indeces
for index_pair in pos_rot_combo:
  pos_index = index_pair[0]; rot_index = index_pair[1]
  #milestone_list.append(raw_milestone_list[index])
  milestone_pos_rot_list.append((raw_milestone_list[pos_index],raw_milestone_list[rot_index]))
  #milestone_pos_rot_list_indeces.append((pos_index,rot_index))

filetree_settings={ # settings for the filetree
  'wet_configs':wet_configs,
  'raw_milestone_list':raw_milestone_list, # because these will need to be updated with the directory information
  'milestone_pos_rot_list':milestone_pos_rot_list,
  #'milestone_pos_rot_list_indeces':milestone_pos_rot_list_indeces,
  #'dry_configs':dry_configs,

}

filetree_settings_all=dict(filetree_settings.items() + sys_params.items())
config_dirlist, md_file_paths, bd_file_paths, raw_milestone_list=filetree.main(filetree_settings_all)
site_list = milestones.split_milestones_by_site(raw_milestone_list)
milestones.write_milestone_file(site_list, milestone_settings['milestone_filename'], master_temperature, md_time_factor, bd_time_factor)

if boolean(inp['ens_equil_colvars']):
  ens_equil_colvars = 'on'
else:
  ens_equil_colvars = 'off'

if sys_params['md']:
  md_settings={ # settings for the md module
    #'milestone_list':milestone_list,
    'milestone_pos_rot_list':milestone_pos_rot_list,
    'raw_milestone_list':raw_milestone_list,
    'configs':wet_configs,
    'receptor_type':inp['receptor_type'],
    'ff':inp['ff'], # the forcefield to use
    'watermodel':inp['watermodel'],
    'cell_shape':inp['cell_shape'],

    'amber_settings': {
      'leap_preload_commands':inp['leap_preload_commands'].split(','),
      'leap_postload_commands':inp['leap_postload_commands'].split(','),
      'leap_program':inp['leap_program'],
      'sample_leap_file':inp['sample_leap_file'],
    },
    'charmm_settings': {
      'parameters':inp['charmm_parameters'].split(','),
      'insert_index':insert_index,
      'recpsf':inp['rec_psf_filename'],
      'ligpsf':inp['lig_psf_filename'],
      'ignore_psf':boolean(inp['ignore_psf'])
    },
    'master_temperature':master_temperature,
    'ligrange':inp['ligrange'],
    'quat_hedron':inp['hedron'],
    'min':boolean(inp['min']), # are we running minimizations at all
    'min_settings':{
      'constrained':inp['min_constrained'], # list what parts of the structure will be constrained during minimizations, including "ligand" (values taken from tcl['lig_indeces'], above), "receptor" (values taken from tcl['rec_indeces']), or a list of all indeces in the pdb file you want constrained
      'restrained':inp['min_restrained'],
      'restrained_force':inp['min_restrained_force'],
      'num_steps':inp['min_num_steps'], # number of minimization steps
      #'out_freq':inp['min_out_freq'],
      'ensemble':inp['min_ensemble'], # irrelevant for minimizations, but necessary for the program
      'namd_settings':{'numsteps':inp['min_num_steps'],'temperature':master_temperature, 'extendedsystem':inp['rec_xsc_filename']},

    },
    'temp_equil':boolean(inp['temp_equil']), # whether we actually run equilibration
    'temp_equil_settings':{
      'constrained':inp['temp_equil_constrained'], # same as above
      'restrained':inp['temp_equil_restrained'],
      'restrained_force':inp['temp_equil_restrained_force'],
      'start_temp':master_temperature,
      'peak_temp':float(inp['temp_equil_peak_temp']),	# These define how the temperature will be adjusted up and then back
      'end_temp':master_temperature,
      'temp_increment':float(inp['temp_equil_temp_increment']),
      'temp_step_time':float(inp['temp_equil_num_steps']), # number of steps per temperature increment
      'ensemble':inp['temp_equil_ensemble'],
      'namd_settings':{'numsteps':inp['temp_equil_num_steps']},
    },
    'ens_equil':boolean(inp['ens_equil']), # whether we will run constrained runs for ensemble equilibrations
    'ensemble_equil_settings':{
      'constrained':inp['ens_equil_constrained'],
      'restrained':inp['ens_equil_restrained'], # whether these structures are harmonically restrained
      'restrained_force':inp['ens_equil_restrained_force'],
      'colvars':boolean(inp['ens_equil_colvars']), # whether collective variables should be imposed between the ligand and the receptor. If set to True, 'constrained' list above must be empty.
      'colvar_settings' : {
        'colvar':inp['ens_equil_colvar_sel'], # list of what parts of the system will have collective variables imposed. Options include 'ligand', 'receptor', 'water', 'relative' (for relative colvars between ligand/receptor), or a list of all indeces in pdb to be constrained
        #'colvar_type':inp['ens_equil_colvar_type'],
        'colvar_force':inp['ens_equil_colvar_force'], # kcal/mol
        'colvarstrajfrequency':inp['ens_equil_colvarstrajfrequency'],
        'colvarsrestartfrequency':inp['ens_equil_colvarsrestartfrequency'],
        'colvar_ligand_indeces':inp['ens_equil_colvar_ligand_indeces'], #map(int, inp['ens_equil_colvar_ligand_indeces'].split(',')),
        'colvar_receptor_indeces':inp['ens_equil_colvar_receptor_indeces'], #map(int, inp['ens_equil_colvar_receptor_indeces'].split(',')),

      },
      'ensemble':inp['ens_equil_ensemble'],
      #'num_steps':'10000', # number of steps to calculate for each ensemble
      'namd_settings':{'numsteps':ens_equil_len,'colvars':ens_equil_colvars,'temperature':master_temperature, 'dcdfreq':dcd_freq,'xstfreq':dcd_freq,'restartfreq':dcd_freq,'outputenergies':dcd_freq, 'outputtiming':dcd_freq, 'veldcdfreq':dcd_freq, 'colvarsconfig':colvars.COLVAR_SCRIPT_NAME},
    },

    'prod_settings':{
      'constrained':inp['fwd_rev_constrained'],

      #'num_steps':'100000',
      'ensemble':inp['fwd_rev_ensemble'],
      'type':inp['fwd_rev_type'], # can be 'protein', 'membrane'
      'namd_settings':{'numsteps':'','tclforces':'on','tclforcesscript':inp['tclforcesscript']},
    },
    'fwd_rev_settings':{
      'constrained':inp['fwd_rev_constrained'],
      'ensemble':inp['fwd_rev_ensemble'],
      'type':inp['fwd_rev_type'], # can be 'protein', 'membrane'
      'namd_settings':{'tclforces':'on','tclforcesscript':inp['tclforcesscript']},
      'dcdfreq':inp['fwd_rev_dcdfreq'],
      'restart_freq':inp['fwd_rev_restart_freq'],
      'run_freq':inp['fwd_rev_run_freq'],
      'xstfreq':inp['fwd_rev_dcdfreq'],
      'extract_xst':inp['extract_xst'],
      'extract_stride':extract_stride, # the number of trajectories to extract from the ensemble simulations
      'extract_first':number_of_ens_equil_frames_skipped, # the number of frames to skip from the beginning of the ensemble trajectory
      'max_num_steps':inp['fwd_rev_max_num_steps'],
      'launches_per_config':inp['fwd_rev_launches_per_config'],
      'frame_chunk_size':inp['fwd_rev_frame_chunk_size']
    },
    'md_file_paths':md_file_paths, # file paths to the MD directories in the anchor file
    'prods_per_anchor':1, # number of simulations per anchor
    #'one_equil_per_anchor':True, # True: all prod simulations will be started from portions of one single equilibration. False: all 50 productions will have their own equilibration
  }


  md_settings['absolute_mode'] = str(absolute_mode)
  md_settings['temp_equil_settings']['temp_range'] = np.concatenate((np.arange(md_settings['temp_equil_settings']['start_temp'],md_settings['temp_equil_settings']['peak_temp'],md_settings['temp_equil_settings']['temp_increment']), \
    np.arange(md_settings['temp_equil_settings']['peak_temp'],md_settings['temp_equil_settings']['end_temp'],-md_settings['temp_equil_settings']['temp_increment']), [md_settings['temp_equil_settings']['end_temp']]))

  md_settings_all = dict(md_settings.items() + sys_params.items() + tcl.items())
  md.main(md_settings_all)

#print "md_settings:", md_settings


if sys_params['bd']:
  
  bd_receptor_dry_pqr=parser.get_structure('bd_receptor_dry_pqr', sys_params['bd_rec_pqr_filename'], pqr=True)
  bd_settings={ # settings for the bd model
    'rec_struct':bd_receptor_dry_pqr,
    'lig_configs': lig_configs,
    'temperature':master_temperature,
    'threads':int(inp['bd_threads']),
    'fhpd_numtraj':inp['fhpd_numtraj'],
    # reaction sites information in milestone_settings below
    'browndye_bin_dir':inp['browndye_bin_dir'],
    'bd_file_paths':bd_file_paths,
    'b_surface_path':os.path.join(sys_params['rootdir'], 'b_surface'),
    'prods_per_anchor':inp['bd_prods_per_anchor'],
    'starting_surfaces':[], # x,y,z,radius
    'ending_surfaces':[], # x,y,z,radius
    'b_surface_ending_surfaces':[], # when the ligand is started from the b-surface, where it can possibly end
    'siteids':[],
    'starting_conditions':'configs', # may be 'spheres' or 'configs'. 'spheres': the ligands are started at random points along the starting_surfaces and then classified into states. 'configs' take configurations from 'lig_configs'
    #'increments':[],
    #'milestones':milestone_list,
    'milestone_pos_rot_list':milestone_pos_rot_list,
    'empty_pqrxml_path':inp['empty_pqrxml_path'],
    'apbs_settings':{
      'apbs_executable':inp['apbs_executable'],
      'ions':[],
      #'ion1rad':inp['ion1rad'], # negative ion
      #'ion2rad':inp['ion2rad'], # positive ion
      #'ion1conc':inp['ion1conc'],
      #'ion2conc':inp['ion2conc'],
      'temp':master_temperature,
    },
    'inputgen_settings':{
      'inputgen_executable':inp['inputgen_executable'],
      'fadd':inp['inputgen_fadd'],
      'gmemceil':inp['inputgen_gmemceil'],
      'resolution':inp['inputgen_resolution'],
      'cfac':inp['inputgen_cfac'],
    },
    
  }
  for key in sorted(inp.keys()):
    if re.match("ion[0-9]+$", key): # then this is an ion
      ion_dict = {'key':key}
      for param in inp[key].split(','): # then pull every string in the list and make a dictionary
        split_param = param.strip().split()
        if len(split_param) < 2: continue
        ion_key = split_param[0]
        ion_value = ' '.join(split_param[1:])
        ion_dict[ion_key] = ion_value
      bd_settings['apbs_settings']['ions'].append(ion_dict)
      
  if inp['ion1conc'] and inp['ion1rad']: # if the old syntax is used
    assert len(bd_settings['apbs_settings']['ions'])==0, "Cannot use 'ion#' parameter at the same time that you are using 'ion1conc' and 'ion1rad' parameters."
    bd_settings['apbs_settings']['ions'].append({'concentration': inp['ion1conc'], 'charge': '-1.0', 'radius': inp['ion1rad'], 'name': 'ion1', 'key': 'ion1'})
    if inp['ion2conc'] and inp['ion2rad']:
      bd_settings['apbs_settings']['ions'].append({'concentration': inp['ion2conc'], 'charge': '1.0', 'radius': inp['ion2rad'], 'name': 'ion2', 'key': 'ion2'})
    
  print "bd_settings['apbs_settings']['ions']:", bd_settings['apbs_settings']['ions']

  #for site in sites: # populate each of the sites for the BD simulations
  for milestone_pair in milestone_pos_rot_list:
    #second_to_last_radius = str(float(site['r']) - float(site['increment']))
    milestone = milestone_pair[0]
    if milestone.bd:
      print "milestone.dimensions:", milestone.dimensions
      bd_settings['starting_surfaces'].append({'x':float(milestone.dimensions['centerx']),'y':float(milestone.dimensions['centery']),'z':float(milestone.dimensions['centerz']), 'outer_radius':float(milestone.dimensions['radius']), 'inner_radius':float(milestone.bd_adjacent.dimensions['radius']), 'inner_index':milestone.bd_adjacent.index, 'outer_index':milestone.index, 'siteid':milestone.siteid})
      #bd_settings['increments'].append(float(milestone.['increment']))
      #bd_settings['ending_surfaces'].append({'x':float(site['x']),'y':float(site['y']),'z':float(site['z']),'radius':float(second_to_last_radius)})
    if milestone.end: # then this is an ending point for our b_surface simulations
      bd_settings['b_surface_ending_surfaces'].append({'x':float(milestone.dimensions['centerx']),'y':float(milestone.dimensions['centery']),'z':float(milestone.dimensions['centerz']),'radius':float(milestone.dimensions['radius']), 'index':milestone.index, 'siteid':milestone.siteid})
    bd_settings['siteids'].append(milestone.siteid)


  #bd_settings = bd_settings.items() + [('milestone_list',milestone_list)])
  bd_settings.update(milestone_settings)

  bd.main(bd_settings)

other_necessary_files={
  'la.tcl':'/path/to/la.tcl',
  'etc':'',
}

# write program paths to a pickle
#print 'namd_special', program_paths['namd_special']
#program_paths_filename=os.path.join(sys_params['rootdir'], 'program_paths.pkl')
#program_paths_file= open(program_paths_filename, 'wb')
#pickle.dump(program_paths, program_paths_file)
#program_paths_file.close

class Test_seekr(unittest.TestCase):
  # several test cases to ensure the functions in this module are working properly
  def test_main(self): # test this function
    print "WARNING: this module does not have comprehensive unit tests"
    return


