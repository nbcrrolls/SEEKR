#!/usr/bin/python

'''
md.py

creates the necessary files to run an MD simulation

'''
import sys, os, math, shutil, subprocess, re #, make_fxd
import transformations
import milestones
import numpy as np
import pdb2 as pdb
from copy import deepcopy # needed to keep track of separate structure objects
import namd_inputs # contains all the details for the namd input files
import colvars
import psfmerge
import  MDAnalysis as mda #import *
from adv_template import Adv_template, File_template
import unittest
import positions_orient
import cPickle as pickle

verbose = True

parser = pdb.Big_PDBParser()

self_path = os.path.dirname(os.path.realpath(__file__)) # get the path to this script

fwd_rev_template_location = os.path.join(self_path, 'fwd_rev.template')

def amber_building(settings):
  '''the pre-minimization procedure for amber ff simulations'''
  i = settings['index']
  building = settings['md_file_paths'][i]['building']
  working_pdb_base = 'holo' #settings['working_pdb_base']
  pdbfile = settings['md_file_paths'][i]['wet_holo']
  amber_settings=settings['amber_settings']
  #sources=amber_settings['sources']
  leap_preload_commands=amber_settings['leap_preload_commands']
  leap_postload_commands=amber_settings['leap_postload_commands']
  leap_program=amber_settings['leap_program']

  leaplist = []
  leapfilename = os.path.join(building,'anchor.leap')
  
  prmtop = os.path.join(building,working_pdb_base+'.parm7')
  inpcrd = os.path.join(building,working_pdb_base+'.rst7')
  newpdb = os.path.join(building,working_pdb_base+'_leap.pdb')

  if os.path.exists(prmtop) and os.path.exists(inpcrd):
    if verbose: print "Amber Parm files already exist. Skipping build phase to save time."
    return prmtop,inpcrd # save us a little time by only running this when it matters

  leapfile = open(leapfilename,'w')

  if amber_settings['sample_leap_file']: # then the user has provided a leap file that we can emulate that should give us what we need
    leaplist = leap_from_sample(settings)
  else:
    for command in leap_preload_commands: # all the commands to be executed before the structure is loaded
      leaplist.append(command.strip())
    leaplist.append('holo = loadpdb %s\n' % (pdbfile,))
    for command in leap_postload_commands: # all the commands to be executed after the structure is loaded
      leaplist.append(command.strip())

    leaplist.append('saveamberparm holo %s %s\n' % (prmtop, inpcrd))
    leaplist.append('savepdb holo  %s\n' % (newpdb,))
    leaplist.append('quit\n')

  # end if block
  leapfile.write('\n'.join(leaplist))
  leapfile.close()


  leapcmd = leap_program+' -f '+leapfilename+' > '+os.path.join(building,'leap.out')
  if verbose: print 'running leap using following command:', leapcmd
  #print 'leapfile', leapfilename
  errcode = 0
  errcode = os.system(leapcmd)
  # check to make sure everything is ok
  if errcode != 0: # then we hit a problem
    errormsg = "LEaP did not run properly. See %s/leap.out for details" % building
    raise Exception, errormsg
  if (not os.path.exists(prmtop)) or (not os.path.exists(inpcrd)):
    errormsg = "LEaP did not generated expected prmtop & inpcrd files. See %s/leap.out for details" % building
    raise Exception, errormsg
  return prmtop, inpcrd

def leap_from_sample(sample_leap_filename, pdbfile, prmtop, inpcrd, newpdb):
  '''Takes a leap file as a sample, along with the other arguments, constructs a list of strings that will be run as a leap program.'''
  sample_file = open(sample_leap_filename, 'r')
  leaplist = []
  print "Now constructing LEAP file from sample file:", pdbfile
  print "If there is an error, make sure that your LEAP file actually works standalone. Then make sure all necessary parameters for all molecules are loaded in the sample LEAP file."
  for line in sample_file.xreadlines():
    line = line.split()
    pdbload_loc = saveamberparm_loc = savepdb_loc = None
    if "loadpdb" in line:
      pdbload_loc = line.index('loadpdb')
    if "saveamberparm" in line:
      saveamberparm_loc = line.index('saveamberparm')
    if "savepdb" in line:
      savepdb_loc = line.index('savepdb')

    if pdbload_loc != None: # then we've found the loadpdb line
      line[pdbload_loc + 1] = pdbfile
    if saveamberparm_loc != None: # then we've found the saveamberparm line
      line[saveamberparm_loc + 2] = prmtop
      line[saveamberparm_loc + 3] = inpcrd
    if savepdb_loc != None: # then we've found the savepdb line
      line[savepdb_loc + 2] = newpdb
    leaplist.append(' '.join(line))
  return leaplist




def charmm_building(settings):
  '''the pre-minimization procedure for charmm ff simulations'''
  # NEED to build new psf file
  if not settings['charmm_settings']['ignore_psf']: # as long as we are not ignoring the creation of psfs
    i = settings['index']
    psf_dir = settings['md_file_paths'][i]['building']
    psf_filename = os.path.join(psf_dir,'holo_wet.psf')
    insert_index = settings['charmm_settings']['insert_index']
    newpsf = psfmerge.merge_psf_files(settings['charmm_settings']['recpsf'], settings['charmm_settings']['ligpsf'], insert_index = insert_index) # merge the psf files
    newpsf.write(psf_filename) # write the new psf files
  return

def parse_selection(selstring):
  """Given a VMD-like selection string, will return a list of atom indeces"""
  try:
    index_list=[]
    raw_list=selstring.strip().split()
    for i in range(len(raw_list)):
      if raw_list[i] == "to": # if using the 'to' keyword, then find the range in between
        index_list += range(int(raw_list[i-1])+1,int(raw_list[i+1])) # get the range from the one previous to the one next
        continue
      else:
        index_list.append(int(raw_list[i]))
  except:
    print "Error: unable to parse index selection text:", selstring
    exit()
  return index_list

def constraints(occ_pdb, constrained, settings): # constraints (atoms in 'constrained' aren't allowed to move
  # Confusingly, constraints in NAMD are different, and should actually be called Restraints. Constraints in NAMD are fixed atoms
  constrained_index_list = []
  if 'ligand' in constrained: # then constrain the ligand
    constrained_index_list += parse_selection(settings['ligrange'])
  if 'receptor' in constrained: # then constrain the receptor
    constrained_index_list += parse_selection(settings['recrange'])
  # if min_constrained contains a list...
  for selection in constrained:
    if type(selection)==type([]): # then it's a list, include
      constrained_index_list += selection

  for atom in occ_pdb.get_atoms():
    if atom.index in constrained_index_list: # changed constrained residue occupancy
      atom.occupancy = 1
    else:
      atom.occupancy = 0

  return

def restraints(occ_pdb, restrained, restrained_force, settings): # harmonic restrains
  # Confusingly, constraints in NAMD are different, and should actually be called Restraints...
  restrained_index_list = []
  if 'ligand' in restrained: # then restrain the ligand
    restrained_index_list += parse_selection(settings['ligrange'])
  if 'receptor' in restrained: # then restrain the receptor
    restrained_index_list += parse_selection(settings['recrange'])
  if 'water' in restrained: # restrain waters
    raise Exception, "restrained waters not yet implemented"
  # if min_constrained contains a list...
  for selection in restrained:
    if type(selection)==type([]): # then it's a list, include all indeces in the list
      restrained_index_list += selection

  for atom in occ_pdb.get_atoms():
    if atom.index in restrained_index_list: # changed restrained residue occupancy
      atom.occupancy = float(restrained_force) # must be converted to a float value
    else:
      atom.occupancy = 0

  return

def make_milestone_list(milestones,hedron):
  '''converts an entire list of milestone dictionaries to TCL-readable array/lists in string form'''
  total_list = ['{']
  for milestone in milestones:

    str_list = ['{']
    if milestone.shape == 'sphere':
      str_list.append('id %i' % milestone.index)
      str_list.append('siteid %s' % milestone.siteid)
      str_list.append('shape sphere')
      str_list.append('center_type %s' % milestone.center_type)
      str_list.append('absolute %s' % milestone.absolute)
      str_list.append('radius %f' % milestone.dimensions['radius'])
      if milestone.center_type == 'coord':
        str_list.append('center { %f %f %f }' % (milestone.dimensions['centerx'], milestone.dimensions['centery'], milestone.dimensions['centerz']))
      elif milestone.center_type == 'atom':
        str_list.append('center_indeces { %s }' % milestone.dimensions['center_indeces'])
        str_list.append('center_com_id ""')
      str_list.append('side outside')
      str_list.append('end %s' % milestone.end)
    elif milestone.shape == 'plane':
      str_list.append('id %i' % milestone.index)
      str_list.append('siteid %s' % milestone.siteid)
      str_list.append('shape plane')
      str_list.append('center_type %s' % milestone.center_type)
      str_list.append('absolute %s' % milestone.absolute)
      str_list.append('distance %f' % milestone.dimensions['distance'])
      str_list.append('normal {%s}' % milestone.dimensions['normal'])
      if milestone.center_type == 'coord':
        str_list.append('center { %f %f %f }' % (milestone.dimensions['centerx'], milestone.dimensions['centery'], milestone.dimensions['centerz']))
      elif milestone.center_type == 'atom':
        str_list.append('center_indeces { %s }' % milestone.dimensions['center_indeces'])
        str_list.append('center_com_id ""')
      str_list.append('side outside')
      str_list.append('end %s' % milestone.end)
    elif milestone.shape == 'rotational':
      str_list.append('id %i' % milestone.index)
      str_list.append('siteid %s' % milestone.siteid)
      str_list.append('shape rotational')
      str_list.append('absolute %s' % milestone.absolute)
      adj_quat_list = make_closest_neighbors(milestone.anchor, hedron)
      str_list.append('quaternion %s' % make_quat_list_in_tcl(milestone.anchor))
      str_list.append( 'adj_quat_list %s' % make_quat_list_in_tcl(adj_quat_list) )
      str_list.append( 'adj_side_list { %s}' % ('unk ' * len(adj_quat_list))) # just make a dummy side list for not, this can be calculated in the TCL
      str_list.append('side uncrossed')
      str_list.append('end %s' % milestone.end)
    str_list.append('}')
    total_list.append(' '.join(str_list))
  total_list.append('}')
  return '\n'.join(total_list)

def make_quat_list_in_tcl(pylist):
  tcllist = []
  if len(pylist) == 0:
    return "{}"
  if type(pylist[0]) == type([]):
    for sublist in pylist:
      tclstr = '{' + ' '.join(map(str,sublist)) + '}'
      tcllist.append(tclstr)
  else:
    tcllist = '{' + ' '.join(map(str,pylist)) + '}'
  return tcllist

def make_closest_neighbors(ref_quat, hedron, tcl_output=True):
  '''given a reference quaternion, will loop through all the quats in hedron, find the ones that tie for closest, and generate a list of the closest quats'''
  # NOTE: this will probably need to be changed once a non-Platonic rotational grid is implemented to account for more distant neighbors
  maxdot = -1.0
  closest_neighbors = []
  for quat in hedron:
    this_dot = np.dot(quat, ref_quat)
    if np.isclose(this_dot, 1.0): continue # skip our own reference quat
    if np.isclose(this_dot, maxdot): # then this might be one of the closest neighbors
      closest_neighbors.append(quat) # then just append this one to the list
    elif this_dot > maxdot: # then our old list doesn't contain the closest neighbors
      closest_neighbors = [quat] # then start the list over
      maxdot = this_dot
  if tcl_output:
    return make_quat_list_in_tcl(closest_neighbors) # convert to TCL before returning
  else:
    return closest_neighbors

def calculate_anticross(cross, anchor):
  '''calculates the opposite rotation than going from anchor to cross'''
  #set CROSS_QUAT $MAIN_NEIGHBOR
  #set anchor_to_cross [quat_mult $CROSS_QUAT [quat_inv $ANCHOR_QUAT]]
  #set anchor_to_anticross [quat_inv $anchor_to_cross]
  #set ANTICROSS_QUAT [ quat_to_anchor_hemisphere  $ANCHOR_QUAT [quat_mult $anchor_to_anticross $ANCHOR_QUAT]]
  anchor_inv = transformations.quaternion_inverse(anchor).tolist() # invert the anchor
  anchor_to_cross = transformations.quaternion_multiply(cross, anchor_inv).tolist() # find quat from anchor to cross
  anchor_to_anticross = transformations.quaternion_inverse(anchor_to_cross).tolist() # invert that quat
  anticross = transformations.quaternion_multiply(anchor_to_anticross, anchor) # get the anticross
  if np.dot(anchor, anticross) < 0.0: anticross *= 1.0 # if its on the wrong hemisphere than the anchor, then flip it back over
  return anticross.tolist()

def close_to_item_in(ouritem, ourlist):
  # sees whether an array has a dot product equal to 1.0 in a list of arrays
  for cmpitem in ourlist:
    if np.isclose(np.dot(ouritem,cmpitem),1.0):
      return True
  return False

def generate_cross_anticross_pairs(ref_quat, allquats):
  '''will find a subset of the closest quats in allquats to refquat'''
  cross_anticross_list = []
  closest_neighbors = make_closest_neighbors(ref_quat, allquats, tcl_output=False)
  cross_list = []
  for neighbor in closest_neighbors:
    anticross = calculate_anticross(cross=neighbor, anchor=ref_quat)
    if not close_to_item_in(anticross, cross_list):
      cross_anticross_list.append(make_quat_list_in_tcl((make_quat_list_in_tcl(neighbor),make_quat_list_in_tcl(anticross))))
      cross_list.append(neighbor)
  return cross_anticross_list

def quat_tclforces(settings):
  '''generates the string that will be included in the NAMD file for the quaternion tclforces script for rotational colvars'''
  tcltemplatestring = """
  set FORCE_CONSTANT {$force_constant}
  set RECROT $recrot
  set REC_INDECES {$rec_indeces}
  set REC_PA1_INDECES {$rec_pa1_indeces}
  set REC_PA3_INDECES {$rec_pa3_indeces}
  set LIG_PA1_INDECES {$lig_pa1_indeces}
  set LIG_PA3_INDECES {$lig_pa3_indeces}
  set LIG_INDECES {$lig_indeces}
  set ANCHOR_QUAT $anchor_quat
  set MAIN_NEIGHBOR $main_neighbor
  set NEIGHBOR_QUATS {
  """
  #for neighbor_quat in settings['neighbor_quats']:
  neighbor_quat_string = ' '.join(settings['neighbor_quats'])
  tcltemplatestring += neighbor_quat_string + '}'
  quat_tcltemplate = Adv_template(tcltemplatestring,settings) # fill in the missing values
  return quat_tcltemplate

def tclforces(settings):
  '''generates the string that will be included in the NAMD file for the tclforces script'''
  tcltemplatestring = """
  set SCRIPT_INTERVAL $script_interval 	;# number of timesteps before script should be evaluated. Examples: 1=script evaluated every timestep. 5=script eval'd every 5th timestep, ...
  set LIGROT $ligrot 					;# whether we care about ligand rotation
  set RECROT $recrot 					;# whether we care about receptor rotation
  set ABORT_ON_CROSSING $abort_on_crossing  ;# whether a job should be aborted right after crossing a milestone
  set PHASE $phase 						;# can be 'forward' or 'reverse'
  set CARE_ABOUT_SELF $care_about_self  ;# whether we care about the milestone's self in this simulation
  set whoami_const $whoami 				;# the index of the milestone we started from
  set whoami $whoami_const
  set whoami_rot_const $whoami_rot 		;# the index of the rotational milestone we started from
  set whoami_rot $whoami_rot_const
  set GRID_EDGE_RAD $grid_edge_rad 		;# radius to the edge of the grid from the center
  set LIGRANGE {$ligrange} 				;# the indeces of the atoms that define the ligand
  set RECRANGE {$recrange} 				;# the indeces of the atoms that define the receptor
  #set RECPA1_LIST {$recpa1_list} 		;# the list of relatively immobile atoms that define the receptor's 1st PA. MARKED FOR REMOVAL
  #set RECPA3_LIST {$recpa3_list} 		;# the list of relatively immobile atoms that define the receptor's 3rd PA. MARKED FOR REMOVAL
  set MILESTONE_LIST $milestone_string  ;# a list of milestones to be converted to arrays
  set SITEID $siteid 					;# the ID if the current site
  set MAX_NUM_STEPS $max_num_steps 		;# the maximum number of steps before starting a new crossing phase
  """
  tcltemplate = Adv_template(tcltemplatestring,settings) # fill in the missing values
  return tcltemplate

def fwd_rev_script(settings):
  script_template = '''
'''

def make_relative_path(oldpath, fromdir):
  '''given a relative 'oldpath' to a file, will construct a relative path from a 'fromdir' directory to the file'''
  oldabspath = os.path.abspath(oldpath)
  fromdirpath = os.path.abspath(fromdir)
  sep = os.path.sep
  oldpathlist = oldabspath.split(sep)
  fromdirlist = fromdirpath.split(sep)
  same = True
  for i in range(len(fromdirlist)): # for every directory on up to the fromdir
    if same and fromdirlist[0]==oldpathlist[0]: # if we are the same up to this point, and these two members are the same
      fromdirlist.pop(0)
      oldpathlist.pop(0)
    elif fromdirlist[0]!=oldpathlist[0]: # as soon as we have a divergence in the path tree
      same=False
      backsteps = len(fromdirlist) # we have to take a few steps back from the 'fromdir'
      backlist = ['..']*backsteps # create a list of '..' backwardses
      relpathlist = backlist + oldpathlist # append the '..'s to the remaining oldpathlist
      break

  return sep.join(relpathlist) # and join them with the '/' to make our relative path

def extract_frames_from_ens(ens_dcd, ens_prmtop, write_dir,start=0, num_frames='all', frame_selection='uniform'):
  '''given the location of a dcd file, will extract some pdb structures from it, write them, and they can be used
  to start other simulations from in the future'''
  # use MDAnalysis to open the file
  universe = mda.Universe(ens_prmtop, ens_dcd)
  allsel = universe.selectAtoms("all") # select all the atoms in here
  traj_len = universe.trajectory.numframes
  if num_frames == 'all': num_frames = traj_len+1-start
  if frame_selection=='uniform':
    frame_list = range(start,(traj_len+1),traj_len/(num_frames-1))[0:num_frames] # we should have a list of length num_frames now, spread uniformly across the
  # elif... # other methods needed, like when it crosses the milestone, which would probably require deeper analysis

  for frame_index in frame_list:
    universe.trajectory[frame_index] # set the frame
    allsel.write(os.path.join(write_dir,"prod%i.pdb" % frame_index)) # write the pdb file to start from
  return frame_list

def get_fhpd(fhpd_file):
  '''given the FHPD filename, returns a list of FHPD trajectory indeces'''
  fhpd_list = []
  for line in open(fhpd_file):
    fhpd_list.append(int(line.strip()))
  return fhpd_list

def prep(settings, holo, stage, inpname, outname='', temperatures=[] ,write_freq='1000',):
  if verbose: print "creating %s files" % stage
  
  i = settings['index']
  path = settings['md_file_paths'][i][stage]
  ff = settings['ff'] # forcefield
  pos_milestone = settings['milestone_pos_rot_list'][i][0]
  rot_milestone = settings['milestone_pos_rot_list'][i][1]
  temperature = settings['master_temperature']
  absolute_mode = settings['absolute_mode']

  if not outname: outname = stage # provide the default outname
  get_cell = False
  if len(temperatures) == 0:
    temperatures = [settings['master_temperature']]
  if stage == 'min':
    stage_settings = settings['min_settings']
    if not stage_settings['namd_settings']['extendedsystem']: # if they didn't define a starting XSC file
      get_cell = settings['cell_shape'] # create the periodic boundary conditions

  if stage == 'temp_equil':
    stage_settings = settings['temp_equil_settings']

  if stage == 'ens_equil':
    stage_settings = settings['ensemble_equil_settings']

  if stage == 'prod':
    stage_settings = settings['prod_settings']

  if stage == 'fwd_rev':
    stage_settings = settings['fwd_rev_settings']

  ensemble = stage_settings['ensemble']
  namd_settings = stage_settings['namd_settings'] # extra settings to send to the NAMD input file
  
  if 'constrained' in stage_settings.keys() and stage_settings['constrained']: # if the settings exists and not empty/False
    fixed = True
    occ_pdb = deepcopy(holo) # copy the holo structure so we can change the occupancy
    constrained = stage_settings['constrained']
    constraints(occ_pdb, constrained, settings) # constrain the selected atoms
    fxd_name = os.path.join(path,'fxd1.pdb')
    namd_settings['fixedAtoms'] = 'on'
    namd_settings['fixedAtomsFile'] = 'fxd1.pdb'
    namd_settings['fixedAtomsCol'] = 'O'
    occ_pdb.save(fxd_name, amber=True, standard=False) # save the constraint file
  else:
    fixed = False

  if 'restrained' in stage_settings.keys() and stage_settings['restrained']: # if the settings exists and not empty/False
    const = True
    occ_pdb = deepcopy(holo) # copy the holo structure so we can change the occupancy
    restrained = stage_settings['restrained']
    restrained_force = stage_settings['restrained_force']
    restraints(occ_pdb, restrained, restrained_force, settings) # constrain the selected atoms
    rest_name = os.path.join(path,'restrained1.pdb')
    occ_pdb.save(rest_name, amber=True, standard=False)
    #ens_namd_settings['fixedatoms'] = 'off' # we may not necessarily want to turn off fixed atoms
    namd_settings['constraints'] = 'on'
    namd_settings['consref'] = 'restrained1.pdb'
    namd_settings['conskfile'] = 'restrained1.pdb'
    namd_settings['conskcol'] = 'O'
  else:
    const = False


  settings['whoami'] = str(pos_milestone.index)
  settings['whoami_rot'] = str(rot_milestone.index)
  settings['siteid'] = str(pos_milestone.siteid)
  settings['grid_edge_rad'] = '0.0' # temporary value until this can be more effectively predicted; will speed up TCL script evaluation
  if ff == 'amber':
    namd_settings['parmfile']=make_relative_path(settings['prmtop'],path) # Find the relative location of the prmtop/inpcrd to be used in the mins
    namd_settings['ambercoor']=make_relative_path(settings['inpcrd'],path)
  elif ff == 'charmm':
    namd_settings['amber'] = 'no'
    namd_settings['coordinates'] = '../holo_wet.pdb'
    namd_settings['structure'] = '../building/holo_wet.psf'
    namd_settings['paratypecharmm'] = 'on'
    counter = 1
    for parameter in settings['charmm_settings']['parameters']:
      if counter == 1: # we need to make sure the namd input file can get more than one parameter file
        i = ""
      else:
        i = str(counter)
      namd_settings['parameters%s' % i] = parameter
      counter += 1
      #namd_settings['parameters2'] = '/extra/moana/rswift1/Permeability/Colvar/A/1F/par_all36_lipid.prm'

  #fhpd_file = 'fhpd.txt'
  namd_settings['watermodel'] = settings['watermodel']
  if stage == 'ens_equil':
    namd_settings['outfilename'] = outname + "_0_1"
    if settings['lig_rot']:
      hedron = positions_orient.get_hedron(settings['quat_hedron'])
      quat_tclforces_settings = {}
      quat_tclforces_settings['force_constant'] = settings['quat_force_constant']
      quat_tclforces_settings['recrot'] = settings['recrot']
      quat_tclforces_settings['rec_indeces'] = settings['quat_rec_indeces']
      quat_tclforces_settings['lig_indeces'] = settings['quat_lig_indeces']
      #quat_tclforces_settings['lig_pa1_indeces'] = settings['ligpa1_list'] # MARKED FOR REMOVAL
      #quat_tclforces_settings['lig_pa3_indeces'] = settings['ligpa3_list'] # MARKED FOR REMOVAL
      #quat_tclforces_settings['rec_pa1_indeces'] = settings['recpa1_list'] # MARKED FOR REMOVAL
      #quat_tclforces_settings['rec_pa3_indeces'] = settings['recpa3_list'] # MARKED FOR REMOVAL
      allquats = make_quat_list_in_tcl(hedron)
      #for i in range(len(allquats)):
      quat = rot_milestone.rotation
      tcl_quat = make_quat_list_in_tcl(quat)
      quat_tclforces_settings['anchor_quat']=tcl_quat
      quat_tclforces_settings['neighbor_quats']=allquats
      closest_neighbors = make_closest_neighbors(ref_quat=quat, hedron=hedron, tcl_output=False)
      if len(closest_neighbors) == 0: closest_neighbors = ["None"]
      number_of_neighbors = len(closest_neighbors)
      #closest_neighbors = make_quat_list_in_tcl(closest_neighbors)
      for f in range(number_of_neighbors):
        quat_tclforces_settings['main_neighbor'] = make_quat_list_in_tcl(closest_neighbors[f])
        namd_settings['tclforces'] = 'on'
        namd_settings['tclforcesscript'] = settings['quatforcesscript'] # the script to be evaluated to keep the ligand forced to the rotational milestone
        namd_settings['tclforces_vars'] = quat_tclforces(quat_tclforces_settings).get_output()
        # write namd file here
        namd_settings['inpfilename'] = inpname
        namd_settings['outfilename'] = outname + "_neighbor_%d_1" % f # 1st number: number of the neighbor, 2nd number: number of ens_equil simulation
        inp, namd_params=namd_inputs.make_input(holo, ff, stage, temperature, write_freq, receptor_type=settings['receptor_type'], ensemble=ensemble, fixed=fixed, get_cell=get_cell, constraints=const, settings=namd_settings) #...
        input_name = os.path.join(path, '%s_%d_1.namd' % (stage,f))
        inp.save(input_name)
        param_filename=(os.path.join(path, 'namd_parameters.pkl'))
        param_file = open(param_filename, 'wb') 
        pickle.dump(namd_params, param_file)
        param_file.close()
      return os.path.join(stage, namd_settings['outfilename']) # return ens_equil
    else: # then lig_rot is False
      namd_settings['inpfilename'] = inpname
      namd_settings['outfilename'] = outname + "_0_1"
      inp, namd_params=namd_inputs.make_input(holo, ff, stage, temperature, write_freq, receptor_type=settings['receptor_type'], ensemble=ensemble, fixed=fixed, get_cell=get_cell, constraints=const, settings=namd_settings) #...
      input_name = os.path.join(path, '%s_%d_1.namd' % (stage,0))
      inp.save(input_name)
      param_filename=(os.path.join(path, 'namd_parameters.pkl'))
      param_file = open(param_filename, 'wb')
      pickle.dump(namd_params, param_file)
      param_file.close()
      return os.path.join(stage, namd_settings['outfilename']) # return ens_equil

  hedron = positions_orient.get_hedron(settings['quat_hedron'])
  allquats = make_quat_list_in_tcl(hedron)
  quat = rot_milestone.rotation
  closest_neighbors = make_closest_neighbors(ref_quat=quat, hedron=hedron, tcl_output=False)
  if len(closest_neighbors) == 0: closest_neighbors = ["None"]
  number_of_neighbors = len(closest_neighbors)

  if stage == 'fwd_rev':
    prelim_string_template = '''set TEMP $temperature
set NUM_REVERSALS $num_reversals
global REV_FILENAME_BASE; set REV_FILENAME_BASE "REV_COMPLETED.txt"
global FWD_FILENAME_BASE; set FWD_FILENAME_BASE "FWD_COMPLETED.txt"
global RESTART_FREQ; set RESTART_FREQ $restart_freq
global RUN_FREQ; set RUN_FREQ $run_freq
set ENS_EQUIL_FIRST_FRAME $begin
set ENS_EQUIL_STRIDE $stride
set LAUNCHES_PER_CONFIG $launches_per_config
set FRAME_CHUNK_SIZE $frame_chunk_size
set UMBRELLA_GLOB_DIR "../ens_equil/" ;# the umbrella sampling directory
set UMBRELLA_GLOB_NAME "ens_equil_0_?.dcd" ;# the umbrella sampling trajs
set nr [numReplicas]
set replica_id [myReplica] ;# get the ID of this replica'''
    hedron = positions_orient.get_hedron(settings['quat_hedron'])
    settings['care_about_self'] = 'False'
    settings['phase'] = 'reverse'
    settings['abort_on_crossing'] = 'True'
    #settings['fhpd_file'] = '../reverse/fhpd.txt'
    settings['max_num_steps'] = stage_settings['max_num_steps']
    settings['milestone_string'] = make_milestone_list(settings['raw_milestone_list'],hedron)
    namd_settings['tclforces_vars'] = tclforces(settings).get_output()

    num_frames = settings['ensemble_equil_settings']['namd_settings']['numsteps']
    dcd_write_freq = int(settings['ensemble_equil_settings']['namd_settings']['dcdfreq'])
    begin = stage_settings['extract_first'] # new parameter
    end = int(num_frames)/dcd_write_freq + 1
    stride = stage_settings['extract_stride'] # new parameter
    launches_per_config = stage_settings['launches_per_config']
    frame_chunk_size = stage_settings['frame_chunk_size']
    restart_freq = stage_settings['restart_freq'] # New parameter
    run_freq = stage_settings['run_freq'] # new parameter
    num_reversals = int((end - begin)/stride)

    namd_settings['inpfilename'] = inpname
    namd_settings['outfilename'] = outname+".$replica_id"
    namd_settings['coordinates'] = "../holo_wet.pdb" # this gets overwritten...
    namd_settings['ambercoor'] = '' # this has to be disabled when coordinates are presented
    namd_settings['velocities'] = ''
    if bool(stage_settings['extract_xst']):
      namd_settings['extendedsystem'] = "../ens_equil/ens_equil_0_1.restart.xsc"
    namd_settings['bincoordinates'] = ""
    namd_settings['binvelocities'] = ""
    namd_settings['temperature'] = temperature
    namd_settings['replicaUniformPatchGrids'] = 'on'
    namd_settings['numsteps'] = ''
    namd_settings['restartfreq'] = '$RESTART_FREQ'
    #namd_settings['id'] = i
    prelim_settings = {'temperature':temperature, 'num_reversals':num_reversals, 'restart_freq':restart_freq, 'run_freq':run_freq, 'begin':begin, 'stride':stride, 'launches_per_config':launches_per_config, 'frame_chunk_size':frame_chunk_size}
    prelim_string = Adv_template(prelim_string_template,prelim_settings).get_output() # fill in the missing values
    namd_settings['prelim_string'] = prelim_string
    post_string_settings = {}
    namd_settings['post_string'] = File_template(fwd_rev_template_location, post_string_settings).get_output()

    inp, namd_params=namd_inputs.make_input(holo, ff, stage, temperature, write_freq, fixed=fixed, ensemble=ensemble, get_cell=get_cell, constraints=const, settings=namd_settings) #...
    input_name = os.path.join(path, 'fwd_rev1.namd')
    inp.save(input_name) 
    param_filename=(os.path.join(path, 'namd_parameters.pkl'))   
    param_file = open(param_filename, 'wb')
    pickle.dump(namd_params, param_file)
    param_file.close()
    #make_fhpd_script(path, 'forward_template.namd', settings['fhpd_file'], number_of_neighbors=number_of_neighbors)
    return os.path.join(stage, namd_settings['outfilename'])


  counter = 1
  for temperature in temperatures:
    namd_settings['inpfilename'] = inpname
    namd_settings['outfilename'] = outname + str(counter)
    inp, namd_params=namd_inputs.make_input(holo, ff, stage, temperature, write_freq, receptor_type=settings['receptor_type'], ensemble=ensemble, fixed=fixed, get_cell=get_cell, constraints=const, settings=namd_settings) #...
    input_name = os.path.join(path, '%s%d.namd' % (stage,counter))
    inp.save(input_name)
    counter += 1
    inpname = namd_settings['outfilename']

  return os.path.join(stage, namd_settings['outfilename'])

def main(settings):
  '''called by seekr, executes other necessary commands'''
  #milestone_list = settings['milestone_list']
  milestone_pos_rot_list = settings['milestone_pos_rot_list']
  for i in range(len(settings['configs'])):
    if milestone_pos_rot_list[i][0].md == False: continue
    settings['index'] = i
    holo = settings['configs'][i]
    ff = settings['ff']
    if ff=="amber":
      if verbose: print "Using forcefield: Amber"
      settings['prmtop'],settings['inpcrd'] = amber_building(settings)
    if ff=="charmm":
      if verbose: print "Using forcefield CHARMM"
      charmm_building(settings)
    last_out=""
    if settings['min']: last_out=prep(settings, holo, stage='min',  inpname='') # if we are supposed to do minimizations
    if settings['temp_equil']: last_out=prep(settings, holo, stage='temp_equil', inpname=os.path.join('..',last_out), temperatures=settings['temp_equil_settings']['temp_range']) # if we are supposed to do temperature equilibrations
    if settings['ens_equil']: last_out=prep(settings, holo, stage='ens_equil', inpname=os.path.join('..',last_out)) # if we are supposed to do temperature equilibrations
    if settings['ensemble_equil_settings']['colvars']: # if running colvars
      colvars.main(settings)

    print "last_out before reverse:", last_out
    #prep(settings,holo, stage='reverse', inpname=os.path.join('..',last_out)) # this can't be run until the ensemble is complete
    #prep(settings,holo, stage='forward', inpname=os.path.join('..',last_out)) # this can't be run until the reverse stage is complete
    prep(settings, holo, stage='fwd_rev', inpname=os.path.join('..',last_out)) # we can put the two stages together



class Test_md(unittest.TestCase):
  # several test cases to ensure the functions in this module are working properly
  def test_main(self): # test this function
    pass

  def test_amber_building(self):
    print "Warning: no unittests written for function 'amber_building' because of it's complicated inputs."
    return

  def test_leap_from_sample(self):
    testfile_name = '/usr/tmp/testfile.leaprc'
    testfile_contents = '''source leaprc.ff03.r1\nsource leaprc.gaff\nloadoff /extra/banzai/lvotapka/projects/seekr/tropc_files/Ca2.lib\nloadoff /extra/banzai/lvotapka/projects/seekr/trypsin_files/benzamidine.lib\nloadamberparams /extra/banzai/lvotapka/projects/seekr/trypsin_files/benzamidine.frcmod\nholo = loadpdb input.pdb\nsaveamberparm holo mine.prmtop mine.inpcrd\nsavepdb holo mine.pdb\ncheck holo\ncharge holo\nquit'''
    testfile = open(testfile_name,'w')
    testfile.writelines(testfile_contents)
    testfile.close()
    output = leap_from_sample(testfile_name, 'test.pdb', 'test.prmtop', 'test.inpcrd', 'testnew.pdb')
    expected = '''source leaprc.ff03.r1\nsource leaprc.gaff\nloadoff /extra/banzai/lvotapka/projects/seekr/tropc_files/Ca2.lib\nloadoff /extra/banzai/lvotapka/projects/seekr/trypsin_files/benzamidine.lib\nloadamberparams /extra/banzai/lvotapka/projects/seekr/trypsin_files/benzamidine.frcmod\nholo = loadpdb test.pdb\nsaveamberparm holo test.prmtop test.inpcrd\nsavepdb holo testnew.pdb\ncheck holo\ncharge holo\nquit'''.split('\n')
    self.assertEqual(output,expected)

  def test_charmm_building(self):
    print "Warning: no unittests written for function 'charmm_building' because of it's complicated inputs."
    return

  def test_parse_selection(self):
    self.assertEqual(parse_selection('1 to 5'), [1,2,3,4,5])
    self.assertEqual(parse_selection('1 3 5'), [1,3,5])
    self.assertEqual(parse_selection('4 to 10 14 16 to 18'), [4, 5, 6, 7, 8, 9, 10, 14, 16, 17, 18])

  def test_constaints(self):
    # using a test pdb file, see if we can properly place constraints
    test_holo = parser.get_structure('small','../test/test_tiny.pdb')
    testsettings = {'ligrange':'1 to 5', 'recrange':'10 to 15'}
    test_consts = "['ligand','receptor']"
    constraints(test_holo, test_consts, testsettings)
    #test_holo.save('/tmp/test_constraints.pdb')
    for i in range(len(test_holo.atoms)):
      if i in range(0,5)+range(9,15):
        self.assertEqual(test_holo.atoms[i].occupancy, 1)
      else:
        self.assertEqual(test_holo.atoms[i].occupancy, 0)
    return

  def test_restraints(self):
    # using a test pdb file, see if we can properly place constraints
    test_holo = parser.get_structure('small','../test/test_tiny.pdb')
    testsettings = {'ligrange':'1 to 5', 'recrange':'10 to 15'}
    test_consts = "['ligand','receptor']"
    restraints(test_holo, test_consts, 3, testsettings)
    #test_holo.save('/tmp/test_constraints.pdb')
    for i in range(len(test_holo.atoms)):
      if i in range(0,5)+range(9,15):
        self.assertEqual(test_holo.atoms[i].occupancy, 3)
      else:
        self.assertEqual(test_holo.atoms[i].occupancy, 0)
    return

  def test_make_milestone_list(self):
    #test_milestone = milestone.Milestones()
    vector = np.array([1,2,3])
    startvector = np.array([0, 0, 0])
    origin = np.array([0,0,0])
    atomid = "2475 2487 2498 2756 2798 2909"
    radius = 2.0; lowest_radius = 2.0; diff_rad = 2.0
    x = 0.0; y = 0.0; z = 0.0
    dimensions = {'center_indeces':atomid, 'radius':radius, 'centerx':x, 'centery':y, 'centerz':z}
    anchor =  tuple(origin + startvector*lowest_radius + vector*diff_rad) # choose the location on the next surface.
    test_milestone = milestones.Milestone(anchor=anchor, dimensions=dimensions, index=0, siteid='site1', absolute=False, center_type='atom')
    expected_result = '{\n{ id 0 siteid site1 shape sphere center_type atom absolute False radius 2.000000 center_indeces { 2475 2487 2498 2756 2798 2909 } center_com_id "" side outside end False }\n}'
    #test_milestones = [{'id':0, }]
    result = make_milestone_list([test_milestone], 'single')
    self.assertEqual(result, expected_result)

  def test_make_quat_list_in_tcl(self):
    testlist = [[1,0,0,0],[0.5,0.5,0.5,0.5],[0,1,0,0]]
    expected_result = ['{1 0 0 0}','{0.5 0.5 0.5 0.5}','{0 1 0 0}']
    self.assertEqual(make_quat_list_in_tcl(testlist),expected_result) # verify equality
    testlist = [1,0,0,0]
    expected_result = '{1 0 0 0}'
    self.assertEqual(make_quat_list_in_tcl(testlist),expected_result) # verify equality

  def test_make_closest_neighbors(self):
    test_ref = [0,1,0,0]
    test_quatlist = [[1,0,0,0],[-0.5,0.5,-0.5,-0.5],[0.5,0.5,0.5,0.5],[0,0,1,0]]
    expected_result = ['{-0.5 0.5 -0.5 -0.5}','{0.5 0.5 0.5 0.5}']
    self.assertEqual(expected_result, make_closest_neighbors(test_ref, test_quatlist))

  def test_calculate_anticross(self):
    print "Warning: no unittests written for function 'calculate_anticross', because it is unused in the current version of SEEKR"

  def test_close_to_item_in(self):
    test_item = [0,1.0,0,0]
    test_list1 = [[1,0,0,0],[-0.5,0.5,-0.5,-0.5], [0,1.0,0,0],[0.5,0.5,0.5,0.5],[0,0,1,0]]
    test_list2 = [[1,0,0,0],[0,0,1.0,0],[-0.5,0.5,-0.5,-0.5],[0.5,0.5,0.5,0.5],[0,0,1,0]]
    result1 = close_to_item_in(test_item, test_list1)
    result2 = close_to_item_in(test_item, test_list2)
    self.assertTrue(result1)
    self.assertFalse(result2)

  def test_generate_cross_anticross_pairs(self):
    print "Warning: no unittests written for function 'generate_cross_anticross_pairs', because it is unused in the current version of SEEKR"

  def test_quat_tclforces(self):
    print "Warning: no unittests written for function 'quat_tclforces', because it is unused in the current version of SEEKR"

  def test_tclforces(self):
    print "Warning: no unittests written for function 'tclforces', due to the complexity of its output"

  def test_make_relative_path(self):
    oldpath = '/usr/bin/gcc'
    fromdir = '/tmp'
    expected_result = '../usr/bin/gcc'
    if not os.path.exists(oldpath) or not os.path.exists(fromdir):
      print "Warning: no unittests run for function 'make_relative_path' due to the fact that the expected paths do not exist on this computer."
    else:
      result = make_relative_path(oldpath, fromdir)
      self.assertEqual(result, expected_result)

  def test_extract_frames_from_ens(self):
    print "Warning: no unittests written for function 'extract_frames_from_ens' because it is unused in the current version of SEEKR"

  def test_get_fhpd(self):
    test_list = [3, 5, 8, 234, 345, 456]
    filename = '/tmp/test_fhpd.txt'
    testfile = open(filename, 'w')
    for number in test_list:
      testfile.write('%s\n' % number)
    testfile.close()
    result = get_fhpd(filename)
    self.assertEqual(test_list, result)

  def test_prep(self):
    print "Warning: no unittests written for function 'prep' because of it's complicated inputs and outputs."

  def test_make_fhpd_script(self):
    print "Warning: no unittests written for function 'make_fhpd_script' because of it's complicated inputs and outputs."

  def test_make_extraction_script(self):
    print "Warning: no unittests written for function 'make_extraction_script' because of it's complicated inputs and outputs."

if __name__ == "__main__":
  #print "Warning: unit tests not yet comprehensive..."
  print "Now testing md.py"
  unittest.main()
