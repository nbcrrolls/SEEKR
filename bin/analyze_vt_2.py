import os, sys, re, random, unittest, glob, tempfile
import argparse
from string import Template
from math import exp, log
from pprint import pprint
import numpy as np
from scipy.sparse import linalg as la
from copy import deepcopy
import xml.etree.ElementTree as ET
from subprocess import check_output
from operator import attrgetter # allows us to sort a list of objects by an object's attribute(s)
import csv
from collections import defaultdict
from itertools import chain
#import pandas as pd
import matplotlib.pyplot as plt
import pickle


k_boltz = 1.3806488e-23
R_GAS_CONSTANT = 0.0019872 # in kcal/mol*K
DEFAULT_TEMP = 0.0 # temperature to assign if not in the xml file
DEFAULT_MD_TIME_FACTOR = 2.0
DEFAULT_BD_TIME_FACTOR = 1000.0
DEFAULT_SHAPE = "sphere"
GRID_TRANSITION = "SEEKR: Milestone Transition: "
GRID_TRANSITION_COMPILE = re.compile(GRID_TRANSITION)
INITIAL_TRANSITION = "source: none,"
INITIAL_TRANSITION_COMPILE = re.compile(INITIAL_TRANSITION)
VT_COLLISION = "SEEKR: Cell Collision: "
VT_COLLISION_COMPILE= re.compile(VT_COLLISION)
INITIAL_COLLISION = "SEEKR: Cell Collision: current: none"
INITIAL_COLLISION_COMPILE = re.compile(INITIAL_COLLISION)

FORWARD_OUTPUT_GLOB = "vt_milestoning_*.out.results"


def boolean(arg):
  if str(arg).lower() in ("false","", "0", "0.0", "[]", "()", "{}"):
    return False
  else:
    return True

class Model():
  ''' object representing the entire milestoning model '''
  def __init__(self): # initialize all variables used in the model
    self.sites = []
    self.temperature = 0.0
    self.md_time_factor = 0.0
    self.bd_time_factor = 0.0
    self.num_sites = 0
    self.num_milestones = 0
    self.positional_milestones = []
    self.rotational_milestones = []
    self.b_surface_milestone = None
    #self.index_dict = {}
    #...

  def add_site(self, site): # append a site to this model
    self.sites.append(site)
    self.num_sites += 1
    #self.num_anchors += site.num_milestones
    #...

  def make_index_dict(self):
    pass

  def make_directories(self): # determines which milestones have data stored in which directories
    i = 0
    for site in self.sites:
      for milestone in site.milestones:
        if milestone.shape == "rotational":
          self.rotational_milestones.append(milestone) # keep track of all rotational milestones for iteration
        else:
          self.positional_milestones.append(milestone)
        i += 1
    if len(self.rotational_milestones) == 0: # if there are no rotational milestones (from, say, an old milestoning xml file)
      self.rotational_milestones.append(milestone(i, 'rotational', False, False, "rotation_%d" % i, anchor="[1, 0, 0, 0]")) # then create one, we need at least one rotational milestone
    for pos_milestone in self.positional_milestones: # for each of the positional milestones
      rot_i = 0
      for rot_milestone in self.rotational_milestones:
        #id_str = '%d_%d' % (site.index, milestone.index)
        directory_name = "anchor_%s_%d" % (pos_milestone.fullname, rot_i)
        print "directory_name:", directory_name

class Site():
  ''' object representing a site: a collection of related milestones that share, say, the same origin '''
  def __init__(self, index, name):
    self.anchors = []
    self.num_anchors = 0
    self.index = index
    self.name = name
    #self.center ?
    #..
  
  def add_anchor(self, anchor):
    self.anchors.append(anchor)
    self.num_anchors +=1

#  def add_milestone(self, milestone):
#    self.milestones.append(milestone)
#    self.num_milestones += 1

class Anchor():
  '''object representing each voronoi cell anchor'''
  def __init__(self, index, md, bd, fullname, directory, siteindex, sitename, coord=None):
    self.milestones = []
    self.num_milestones = 0
    self.index = index
    self.fullname = fullname
    self.directory = directory
    self.coord = coord
    self.bd = bd
    self.md = md
    self.site = siteindex
    self.sitename = sitename
    self.total_steps = 0
    self.offsets = {}

  def add_milestone(self, milestone):
    self.milestones.append(milestone)
    self.num_milestones += 1

  def parse_md_transitions(self, ):
    'find all forward phase simulation output files'
    forward_dir_glob = os.path.join(self.directory,'md','fwd_rev',FORWARD_OUTPUT_GLOB)
    #print forward_dir_glob
    forward_output_filenames = sorted(glob.glob(forward_dir_glob))
    # read files and sift out the transition lines
    unsorted_transitions = []
    unsorted_collisions = []
    replicate = 1
    self.offsets[replicate] = 0
    #print 'current max step' , info['max_steps']
    for filename in forward_output_filenames:
      transition_lines = []
      vt_collisions = []
      file_max_steps = 0  
      print 'parsing transitions from file:', filename
      for line in open(filename,'r'):
        if re.match(GRID_TRANSITION_COMPILE, line):
          if re.search(INITIAL_TRANSITION_COMPILE, line): continue #we don't want to count a transition for the first milestone touched by the simulation  
          else:
            transition_lines.append(line)
        if re.match(VT_COLLISION_COMPILE, line):
          if re.match(INITIAL_COLLISION_COMPILE, line): continue
          else:
            vt_collisions.append(line)
           # print line
      #pprint(transition_lines)
    # feed the lines into the Transition object
      for line in transition_lines:
        transition = Transition(line)
        transition.replicate = replicate
        unsorted_transitions.append(transition)

      for line in vt_collisions:
        collision = Collision(line)
        collision.replicate = replicate
        unsorted_collisions.append(collision)
      file_max_steps = collision.step
      
      self.total_steps += file_max_steps
      self.offsets[replicate+1] = self.offsets[replicate]+file_max_steps 
      replicate += 1
    self.transitions = unsorted_transitions
    self.collisions = unsorted_collisions
    print 'anchor', self.index, self.offsets
    return self.total_steps

  def get_md_transition_statistics(self, md_time_factor=DEFAULT_MD_TIME_FACTOR, max_step=None, ):
    'parse the transition data to obtain transition counts'
    counts = {} # all the sources and their destinations
    #cell_counts = {} # keeps track of transitions between voronoi cells
    #cell_times = {} #keeps track of total time simulated ina voronoi cell
    total_counts = {} # keeps track of all counts to any destination
    total_times = {} # the total time of all counts, to be averaged later
    avg_times = {} # the times to transition out of each source
    site_index = self.site
    site_name = self.sitename
   
    for transition in self.transitions:
 
      source = transition.src
      dest = transition.dest
      time = transition.time
      anchor = transition.anchor
      stepnum = transition.cur_step
      src_key = source
      dest_key = dest

      if max_step != None and int(stepnum + self.offsets[transition.replicate]) > max_step:
        break
      if self.index in counts.keys():
        if src_key in counts[self.index].keys():
          if dest_key in counts[self.index][src_key].keys():
            counts[self.index][src_key][dest_key] += 1
          else:
            counts[self.index][src_key][dest_key] = 1
          total_times[self.index][src_key] += time * md_time_factor
          total_counts[src_key] += 1
        else:
          counts[self.index][src_key] = {dest_key:1}
          total_counts[src_key] = 1
          total_times[self.index][src_key] = (time * md_time_factor)
      else:
        counts[self.index]= {src_key: {dest_key: 1}}
        total_times[self.index] = {src_key : (time * md_time_factor)}
        total_counts[src_key] =1

    return counts, total_counts, total_times, avg_times

  def get_md_vt_collisions(self, md_time_factor=DEFAULT_MD_TIME_FACTOR, max_step=None,):
    
    cell_counts = {}
    for collision in self.collisions:
      if max_step != None and int(collision.step + self.offsets[collision.replicate]) > max_step:
        break

      curr_cell = collision.curr_cell
      new_cell = collision.new_cell   
      if curr_cell in cell_counts.keys():
        if new_cell in cell_counts[curr_cell].keys():
          cell_counts[curr_cell][new_cell] += 1
        else:
          cell_counts[curr_cell][new_cell] = 1
      else:
        cell_counts[curr_cell] = {new_cell:1} 
    total_time = collision.step *md_time_factor + self.offsets[collision.replicate] * md_time_factor 

    return cell_counts, total_time

class Milestone():
  ''' object representing a single milestone '''
  def __init__(self, id, shape, end, normal=None, radius=None, ):
    self.id = id # most of this information is read from the provided milestone.xml file
    self.shape = shape
    self.end = end.lower()
    #self.fullname = fullname
    #self.directory = directory
    #self.anchor = anchor
    self.normal = normal
    self.radius = radius
    #self.bd = bd
    #self.md = md
    #self.site = siteindex
    #self.sitename = sitename
    self.transitions = [] # all the transition statistics associated with this milestone

class Collision():
  '''Object representing all collisiont with a voronoi cell edge'''
  def __init__(self,line):
    line = line.strip() # we have to parse the line    
    self.line = line
    linetail = line[len(VT_COLLISION)-1:] # just take the last bit of the line, the important part, but not the endline
    linelist = linetail.split(',') # split line into a list of elements
    dictlist = map(lambda a: a.strip().split(': '), linelist) # map the line to a list of lists for dictionary conversion
    linedict = dict(dictlist) # convert the list of lists into a dictionary
    #self.src = int(linedict['current'].strip())
    self.curr_cell = int(linedict['current'].strip())
    self.new_cell = int(linedict['new'].strip())
    self.step = int(linedict['stepnum'].strip())

class Transition():
  ''' object representing a transition event between one milestone to another. Used to construct transition statistics'''
  def __init__(self, line):
    line = line.strip() # we have to parse the line
    self.line = line
    linetail = line[len(GRID_TRANSITION)-1:] # just take the last bit of the line, the important part, but not the endline
    linelist = linetail.split(',') # split line into a list of elements
    dictlist = map(lambda a: a.strip().split(': '), linelist) # map the line to a list of lists for dictionary conversion
    linedict = dict(dictlist) # convert the list of lists into a dictionary
    self.src = int(linedict['source'].strip())
    self.dest = int(linedict['destination'].strip())
    self.cur_step = float(linedict['stepnum'].strip())
    self.time = float(linedict['incubation time'].strip().split()[0])
    self.anchor = int(linedict['anchor'].strip())
    #self.curr_cell = int(linedict['curr_cell'].strip())
    #self.new_cell = int(linedict['new_cell'].strip())
    #self.ligand_com = linedict['ligand COM'].strip().split()
    #self.receptor_com = linedict['receptor COM'].strip().split()
    #self.receptor_start_com = linedict['receptor start COM'].strip().split()
    #self.umbrella_step = int(linedict['ID'].strip().split()[0])
    #self.velocity_step = int(linedict['VEL_ID'].strip().split()[0])

  def print_status(self):
    print "src:", self.src
    print "dest:", self.dest
    print "cur_step:", self.cur_step
    print "time:", self.time
    #print "ligand_com:", self.ligand_com
    #print "receptor_com:", self.receptor_com
    #print "receptor_start_com:", self.receptor_start_com

def parse_milestoning_file(milestoning_filename):
  'given a milestoning file, will parse the XML and generate a model object'
  tree = ET.parse(milestoning_filename) # now the entire XML file is in a dictionary
  root = tree.getroot() # object at bottom of the tree
  model = Model() # create the model object

  # Temperature
  xml_temp = root.find('temperature') # temperature tag
  if xml_temp != None: # make sure it exists
    model.temperature = float(xml_temp.text.strip())
  else: # some old milestones.xml files may not have this tag
    model.temperature = DEFAULT_TEMP

  # MD Time Factor
  xml_md_time_factor = root.find('md_time_factor') # temperature tag
  if xml_md_time_factor != None: # make sure it exists
    model.md_time_factor = float(xml_md_time_factor.text.strip())
  else: # some old milestones.xml files may not have this tag
    model.xml_md_time_factor = DEFAULT_MD_TIME_FACTOR

  # MD Time Factor
  xml_bd_time_factor = root.find('bd_time_factor') # temperature tag
  if xml_bd_time_factor != None: # make sure it exists
    model.bd_time_factor = xml_bd_time_factor.text.strip()
  else: # some old milestones.xml files may not have this tag
    model.xml_bd_time_factor = DEFAULT_BD_TIME_FACTOR


  site_counter = 0 
  for branch in root:
    if branch.tag != "site": continue # make sure that we are parsing a site
    site_name = branch.find('name').text.strip()
    print "Now parsing milestones from site:", site_name, "in XML file."
    site_obj = Site(site_counter, site_name) # create the site object
    for anchor in branch: #iterate through each anchor in the site
      if anchor.tag != "anchor": continue #ensure we are reading voronoi anchors
      index = anchor.find('index').text.strip()
      coord = anchor.find('coord').text.strip()
      fullname = anchor.find('fullname')
      if fullname != None:
        fullname = fullname.text.strip() # parse all the attributes of the milestone
      else: # then just name it by its anchor
        fullname = str(anchor)
      directory_text = anchor.find('directory').text
      if directory_text:
        directory = directory_text.strip() # directory where this milestone is located
      else:
        directory = None
      bd =boolean(anchor.find('bd').text.strip())
      md =boolean(anchor.find('md').text.strip())
      anchor_obj = Anchor(index, md, bd, fullname, directory, site_counter, site_name, coord) 
        

      for milestone in anchor: #read all of the milestones in this anchor
        if milestone.tag != "milestone": continue
        id = milestone.find('id').text
        radius = None
        normal = None
        shape_xml = milestone.find('shape')
        if shape_xml != None:
          shape = shape_xml.text.strip()
        else:
          shape = DEFAULT_SHAPE # by default
        if shape == "plane": # get other information based on the milestone shape
        #  anchor = milestone.find('anchor').text.strip()
          normal = milestone.find('normal').text.strip()
        elif shape == "sphere":
        #  anchor = milestone.find('anchor').text.strip()
          radius = milestone.find('radius').text.strip()
        elif shape == "rotational":
          pass
        end = milestone.find('end').text.strip()
        milestone_obj = Milestone(id, shape, end, normal, radius)
        anchor_obj.add_milestone(milestone_obj)
      site_obj.add_anchor(anchor_obj)
      
    model.add_site(site_obj)
    site_counter += 1 

      #model.b_surface_milestone = Anchor(index="0", shape="sphere", end="True", md=False, bd=True, fullname="b_surface", directory="b_surface", siteindex=0, sitename="b_surface")

  return model 

def add_dictionaries(dict1, dict2):
  '''
  adds the values numerically within each dictionary
  NOTE: dict1 is updated and returned BY REFERENCE
  '''
  new_dict = dict1
  for key in dict2.keys():
    if key in dict1.keys():
      dict1[key] += dict2[key]
    else:
      dict1[key] = dict2[key]

  return dict1

def parse_bd_results(bd_results_filename):
  ''' given a BD results file name, will open the file and extract information about state transitions'''
  #bd_results_file = open(bd_results_filename, 'r')
  bd_dict = {}
  bd_time = None
  tree = ET.parse(bd_results_filename)
  root = tree.getroot()
  for tag in root:
    if tag.tag == "reactions":
      reactions = tag
      for tag2 in reactions:
        i = 0
        if tag2.tag == "escaped":
          bd_dict['inf'] = int(tag2.text)
        elif tag2.tag == "completed":
          site = tag2[0].text.strip()
          #print "site:", site
          n = tag2[1].text
          #print "n:", n
          #name = outer_state[i] + '_' + str(site)
          bd_dict[site] = int(n)
          i += 1
        elif tag2.tag == "time":
          bd_time = float(tag2.text)

  return bd_dict, bd_time

def parse_bound_state_args(bound_args):
  bound_dict = {}
  bound_pairs = bound_args.split(',')
  for pair in bound_pairs:
   # print 'PAIR'+ pair
    site_index = pair.split(':')
   # print site_index
    if len(site_index) == 1:
      site = 'all'
      index = site_index[0]
    elif len(site_index) == 2:
      site = site_index[0]
      index = site_index[1]
    if site not in bound_dict:
      bound_dict[site] = [index]
    else:
      bound_dict[site].append(index)
   # print bound_dict
  return bound_dict


def read_transition_statistics_from_files(model, verbose):
  '''This function parses the transitions statistics from the simulation output files for later analysis'''
  #info = {'max_steps':0, }
  #info = {}
  total_steps = 0
  for site in model.sites:
    for anchor in site.anchors:
    #for milestone in site.milestones:
      if anchor.md == True and anchor.directory:
        #if verbose: print 'parsing md transitions for:Anchor', milestone.fullname
        #print info['max_steps']
        print 'parsing md transitions for:Anchor', anchor.fullname
        max_steps = anchor.parse_md_transitions()
        print max_steps, total_steps
        if max_steps > total_steps:
          total_steps = max_steps
        #else:
        #  print "last anchor longer"
  
  return total_steps 

def analyze_kinetics(calc_type, model, bound_dict, doing_error, verbose, bd_time, md_time_factor = DEFAULT_MD_TIME_FACTOR, max_steps =None, error_number = 1, error_skip = 1, milestone_conv='False' ,  conv_stride=50, rand_conv= 'False', rand_samples=10, plt_name= ''):
  '''main function to perform all kinetics analyses.
  Given a Model() object with its statistics filled out, it will return an estimate of the kinetic
  value, given all or a subset of its statistics.
  '''
  counts = {}; times = {}; total_counts = {}; total_cell_counts = {}; total_times = {}; avg_times = {}; trans = {}; total_cell_times = {}; T_a = {}
  end_indeces = [];
  N = {}
  if verbose: print 'max_steps', max_steps 
  for site in model.sites:
    for anchor in site.anchors:
      if anchor.md == True and anchor.directory:
        if verbose: print 'Anchor', anchor.fullname
        this_counts, this_total_counts, this_total_times, this_avg_times = anchor.get_md_transition_statistics(model.md_time_factor, max_steps)
        this_cell_counts, this_cell_time = anchor.get_md_vt_collisions(model.md_time_factor, max_steps)
        total_counts = add_dictionaries(total_counts, this_total_counts)
        if verbose: print 'counts',  this_counts
        total_cell_counts = add_dictionaries(total_cell_counts, this_cell_counts)
        total_cell_times[int(anchor.index)] = this_cell_time
        if verbose: print 'times', this_total_times
        if verbose: print 'cell times', total_cell_times
        #total_times = add_dictionaries(total_times, this_total_times)
        #total_cell_times = add_dictionaries(total_cell_times, this_cell_times)
        for src_key in this_counts.keys():
          if src_key in counts.keys():
            counts[src_key] = add_dictionaries(counts[src_key], this_counts[src_key])
          else:
            counts[src_key] = this_counts[src_key]
        #print "len(transitions)", len(milestone.transitions)
      for milestone in anchor.milestones:
        if milestone.end == "true": # then its an end milestone, and there will be no transitions out of it
          end_indeces.append(int(milestone.id))
      for src_key in this_total_times.keys():
        if src_key in times.keys():
          times[src_key] = add_dictionaries(times[src_key], this_total_times[src_key])
        else:
          times[src_key] = this_total_times[src_key]


#      if anchor.bd == True and anchor.directory:
#        this_counts, this_total_counts, this_total_times, this_avg_times = milestone.get_bd_transition_statistics(bd_time=bd_time)
#        print 'TIME', this_avg_times
#        total_counts = add_dictionaries(total_counts, this_total_counts)
#        total_times = add_dictionaries(total_times, this_total_times)
#        for src_key in this_counts.keys():
#          if src_key in counts.keys():
#            counts[src_key] = add_dictionaries(counts[src_key], this_counts[src_key])
#          else:
#            counts[src_key] = this_counts[src_key]

 # for src_key in total_times.keys(): # construct the average incubation times
 #   R_cell[src_key] = total_times[src_key] / total_cell_times[src_key]

#  for src_key in counts.keys():
#    temp = {}
#    for dest_key in counts[src_key].keys(): # make the transition probability dictionary
      #temp[dest_key] = float(counts[src_key][dest_key]) / float(total_counts[src_key])
      #N_cell[src_key][dest_key]



## Calculate Voronoi cell equilibrium probability ##
  #for cell_src_key in total_cell_counts.keys():
  #total_cell_times = {0: 84668000, 1: 86494000, 2:86504000, 3:92662000, 4:43832000, 5:42430000, 6:43522000} #hard coded -- includes times from combined face sampling
  #print "total_cell_counts"; pprint(total_cell_counts)  
  #print "total times"; pprint(total_times)
  #print "cell times"; pprint(total_cell_times)
  k_cell = np.zeros((len(total_cell_times),len(total_cell_times)))  
  k_mod = np.zeros((len(total_cell_times),len(total_cell_times))) 
  



  for cell in total_cell_counts.keys():
    for new_cell in total_cell_counts[cell].keys():
      if new_cell == -1: continue #skip transitions to bound state milestone
      elif new_cell not in end_indeces: #hard code for testing
        #print cell, ", ",  new_cell
        #print total_cell_counts[cell][new_cell]
        #print total_cell_times[cell]
        k_cell[cell][new_cell] = (float(total_cell_counts[cell][new_cell])/float(total_cell_times[cell]))
        #print k_cell[cell][new_cell]
     
  #print "k_cell"; pprint(k_cell)
 
## Create the steady state flux matrix##

  for i in range(len(k_cell[1])):
    for j in range(len(k_cell[1])):
      if i == j:
        k_mod[i][j] = -(np.sum(k_cell[i]))
      else:
        k_mod[i][j] = k_cell[j][i]

## Substitute redundant equation with normalization condition

  k_mod[-1][:] = 1
  #print "k_mod"; pprint(k_mod)

  p_0= np.zeros((len(total_cell_times)), dtype="float")  
  p_0[-1] = 1.0

  #pprint(p_0)

  #print "k_cell:", np.shape(k_cell)
  #print "p_0:", np.shape(p_0)
  #k_cell_trans = np.transpose(k_cell)




## Calculate the equilibrium probabilities for each voronoi cell
  p_equil = np.linalg.solve(k_mod, p_0)
  if verbose: pprint(p_equil)
  #print np.shape(p_equil)
  

  p_equil_ref = p_equil[-1]
  #print p_equil[-1]
  #print p_equil_ref
  #print range(len(p_equil))
  delta_G = np.zeros(len(p_equil))
  for i in range(len(p_equil)):
    #print i
    #print p_equil[i] 
    delta_G[i] = -model.temperature * R_GAS_CONSTANT * log(p_equil[i] / p_equil_ref)
    #print delta_G[i]


  print "Delta G: "; pprint(delta_G)
## Using the V cell equilibrium probabilities, calculate the rate matrix, Q
  if verbose: print "counts: ", counts
  if verbose: print "times: ", times

  T_a = np.zeros(len(p_equil))
  for cell in total_cell_times:
    T_a[cell] = p_equil[cell]/total_cell_times[cell]
  T_tot = 1/np.sum(T_a)


  N = np.zeros((len(p_equil)+1,len(p_equil)+1))
  N_conv = np.zeros((len(p_equil)+1,len(p_equil)+1))
  for anchor in counts.keys():
    for src in counts[anchor].keys():
      for dest in counts[anchor][src].keys():
        #print p_equil[int(anchor)]
        #print counts[anchor][src][dest]
        #print total_cell_times[int(anchor)]
        N[src][dest] = p_equil[int(anchor)] * float(counts[anchor][src][dest])/ total_cell_times[int(anchor)]
        if max_steps != None: 
          if  total_cell_times[int(anchor)] >= max_steps * md_time_factor: 
            N_conv[src][dest] = float(counts[anchor][src][dest])/ total_cell_times[int(anchor)]
          else:
            N_conv[src][dest] = np.nan      

  if verbose: print "N:", N

  R = np.zeros(len(p_equil)+1)
  R_conv = np.zeros((len(p_equil)+1,len(p_equil)+1))
  for anchor in times.keys():
    for src in times[anchor].keys(): 
      R[src] += (p_equil[int(anchor)] * times[anchor][src]/ total_cell_times[int(anchor)])
      if max_steps != None:
        if total_cell_times[int(anchor)] >= max_steps* md_time_factor:
          R_conv[int(anchor)][src] = times[anchor][src]/ total_cell_times[int(anchor)] 
        else:
          R_conv[int(anchor)][src] = np.nan

  if verbose: print "R:", R


  Q = np.zeros((len(p_equil)+1,len(p_equil)+1))
  for i in range(len(N[0])):
    for j in range(len(N[0])):
      Q[i][j] = N[i][j]/R[i]
      

  for i in range(len(N[0])):
    Q[i][i] = -np.sum(Q[i])



  if verbose: print ""
  if verbose: print Q

# Calculate MFPT and off rate
  
  T= calc_MFPT_vec(Q)

  #MFPT = T[0]
  #k_off = 1e15/MFPT
  #print "k_off = ", k_off, "s-1" 

  total_sim_time = 0
  for i in total_cell_times.keys():
    print i, total_cell_times[i]/1e6, "ns"
    total_sim_time += total_cell_times[i]

  print "Total simulation time: " ,  total_sim_time/1e6, "ns" 

  return p_equil, N, R, T, T_tot, Q, N_conv, R_conv, k_cell, 

def calc_MFPT_vec(Q):
  Q_hat = Q[:-1,:-1]

  #if verbose: print Q_hat

  I = np.zeros(len(Q_hat[0]))
  #I[:] = np.sqrt(len(Q_hat[0]))
  I[:] = 1
  
  T= np.linalg.solve(Q_hat,-I)

  #MFPT = T[0]
  #if verbose: print "MFPT =", T, "fs"

  #k_off = 1e15/MFPT

  #print "T", T  

  return T

def monte_carlo_milestoning_error(Q0, N, R, p_equil, T_tot, num = 1000, skip = 0):
  '''Samples distribution of rate matrices assumming a poisson (gamma) distribution with parameters Nij and Ri using Markov chain Monte Carlo
    Enforces detailed Balance-- using a modified version of Algorithm 4 form Noe 2008 for rate matrices.--  
  Distribution is:  p(Q|N) = p(Q)p(N|Q)/p(N) = p(Q) PI(q_ij**N_ij * exp(-q_ij * Ri))
  '''
  m = N.shape[0] #get size of count matrix
  Q = Q0
  Q_mats = []

  print "Q", Q.shape
  print Q
  P = np.zeros((m,m))
  Q_test = np.zeros((m,m))
  tau = np.zeros(m)

## test for convirting Q to P and tau
  for i in range(m):
    for j in range(m):
      if i == j: continue 
      Q_test[i,j] = Q[i,j]


  #Q_test[:][-1] = 1 
  #print "Q_test", Q_test.shape, Q_test
  #print "transpose", Q_test.T

  
  # foo = np.zeros(m, dtype="float")
  # foo[-1] = 1.0
  # foo2 = np.transpose(foo)

  # print "foo", foo2.shape




  for i in range(m):
    for j in range(m):
       if i==j: continue
       P[i,j] = Q_test[i,j] / np.sum(Q_test[i])
    #tau[i] = 1 / np.sum(Q_test[i])
  
  print "P" , P
  # print "tau", tau

  #print "p_equil"
  #print p_equil

  ## calculate initial equilibrium flux by solving: pi = pi Q0  -- left eigenvector of Q0
  val, vec = la.eigs(P.T, k = 1, which= 'LM')

  #print val.real
  #print "vec", vec.real

  pi = vec[:,0].real

  #print "pi", pi
  pi /= pi.sum()

  #print "normalized", pi

  #print sum(pi)

  pi_ref = pi[-1]
  dg = np.zeros(m)


  #for i in range(m):
  #  dg[i] = -300 * R_GAS_CONSTANT * log(pi[i] / pi_ref)

  #print "milestone Free energies", dg
  #print "new_vec", new_vec
  #pi = np.real(vec[:, idx])
  # pi/pi.sum()

  # P_hat = P[:-1,:-1]
  # print "p hat", P_hat

  # foo = np.zeros(m-1, dtype="float")
  # foo[:] = 1.0
  # #foo = np.transpose(foo)
  # P_mod = np.zeros((m-1,m-1))
  # for i in range(m-1):
  #   for j in range(m-1):
  #     P_mod[i,j] = P_hat[i,j]
  # P_mod[-1][:] = 1.0  #sub redundant equation with normalization condition
  # print "p mod", P_mod.shape, P_mod
  # print foo.shape, foo


  # pi= np.linalg.solve(P_mod,foo)

  # print "pi", pi

  # print "foo"

  for counter in range(num*(skip+1)):

    ##Generate random variables for step
    r1 = random.random() #genarate uniform random number for do reversible element shift
    r2 = random.random() #genarate uniform random number for acceptance probability
    Qnew = Q
    print "r1", r1

    if r1 < 0.5:  #then do reversible element shift
      i = random.randint(0,m-1)
      j = random.randint(0,m-1)
      print "r" , r2, "i", i, "j", j 
      #print "N", N.shape, N
      #print "R", R.shape, R
      #print "Q new", Qnew


      if Qnew[i][j] == 0.0: 
        print "skip"
        continue
      if i == j: 
        print "skip"
        continue
      

      print max(Qnew[i,i], pi[j]/pi[i] * Qnew[j,j]), Qnew[i,j]
      ## Generate random perturbation, delta, bounded to prevent any Q from changing sign
      #delta = random.uniform(max(Qnew[i,i], (pi[j]/pi[i] * Qnew[j,j])), Qnew[i,j])
      delta = random.uniform(-1E99 , min(-Qnew[i,i], -(pi[j]/pi[i] * Qnew[j,j]), Qnew[i,j]))
      print "delta", delta

      print np.sqrt((Qnew[i,j] - delta)**2 + (Qnew[j,j] - (pi[i]/pi[j]) * delta)**2 /(Qnew[i,j]**2 +Qnew[j,i]**2))
      print ((Qnew[i,j] - delta)/Qnew[i,j])**N[i,j]
      print ((Qnew[j,i] - (delta * pi[i]/pi[j])) / Qnew[j,i])**N[j,i]
      print Qnew[j,i], pi[i], pi[j], N[j,i]
      print pi[i]/pi[j]
      print delta * (pi[i]/pi[j])
      print (Qnew[j,i] - (delta * pi[i]/pi[j])) / Qnew[j,i]


      #print ((Qnew[i,i] + delta) / Qnew[i,i])**N[i,i]
      #print ((Qnew[j,j] + (pi[i]/pi[j]) * delta) / Qnew[j,j])**N[j,j]

      print "pi", pi

      #Calculate acceptance probability for reversible element shift
      p_acc = np.sqrt((Qnew[i,j] - delta)**2 + (Qnew[j,j] - (pi[i]/pi[j]) * delta)**2 /(Qnew[i,j]**2 +Qnew[j,i]**2)) * ((Qnew[i,j] - delta)/Qnew[i,j])**N[i,j] * ((Qnew[j,i] - (pi[i]/pi[j] * delta)) / Qnew[j,i])**N[j,i] #* ((Qnew[i,i] + delta) / Qnew[i,i])**N[i,i] * ((Qnew[j,j] + (pi[i]/pi[j]) * delta) / Qnew[j,j])**N[j,j]
      print "p_acc", p_acc
      
      if r2 <= p_acc:
        print "performing reversible element shift..."
        

        Qnew[i,i] = (Qnew[i,i]) + delta
        Qnew[j,j] = (Qnew[j,j]) + (pi[i] / pi[j] * delta)
        Qnew[i,j] = Qnew[i,j] - delta
        Qnew[j,i] = Qnew[j,i] - (pi[i] / pi[j] * delta)
        # Qnew[i,i] = -(abs(Qnew[i,i]) + delta)
        # Qnew[j,j] = -(abs(Qnew[j,j]) + (pi[i] / pi[j] * delta))
        # Qnew[i,j] = abs(Qnew[i,j]) - delta
        # Qnew[j,i] = abs(Qnew[j,i]) - (pi[i] / pi[j] * delta)

        print Qnew

      # else: # do row shift
      #   #generate random variables
      #   print "check for row shift..."
      #   i = random.randint(0,m-1)
      #   alp = random.uniform(0,(1/(1-Q[i,i])))
      #   print "alpha", alp

      #   #print np.sum(N[i])

      #   p_acc = alp ** (m - 2 + np.sum(N[i]) - N[i,i]) * ((1 - alp * (1-Qnew[i,i]))/Qnew[i,i]) ** N[i,i]
    

      #   print "p_acc 2", p_acc

      #   if r2 <= p_acc:
      #     print "performing row shift"
      #     for j in range(0,m-1):
      #       Qnew[i,j] = alp * Qnew[i,j]
      #     Qnew[i,i] = 0
      #     Qnew[i,i] = 1- np.sum(Qnew[i])

      #     print Qnew

      #     ## Update stationary distribution
      #     print "updating stationary distribution"

      #     for j in range(m):
      #       pi[j] = alp * pi[j] /pi[i] + alp *(1-pi[i])

      #     pi[i] = 0
      #     pi[i] = 1- np.sum(pi)

      #     print pi

    #print "Pi", pi 
    if counter % (skip +1) == 0:
      Q_mats.append(Qnew)
    Q = Qnew
  return Q_mats


def monte_carlo_milestoning_nonreversible_error(Q0, N_pre, R_pre, p_equil, T_tot, num = 20, skip = 0):
  ''' Samples a distribution of rate matrices that are nonreversible
      using Markov Chain Monte Carlo method

      The distribution being sampled is:
      p(Q|N) = p(Q)p(N|Q)/p(N) = p(Q) PI(q_ij**N_ij * exp(-q_ij * R_ij))

      N = count matrix
      R = transition times

  '''
  #Q0, N_sum = count_mat_to_rate_mat(N, avg_t) # get a rate matrix and a sum of counts vector
  m = Q0.shape[0] # the size of the matrix
  Q = Q0
  Q_mats = []
  N = []
  R = []
  
  N = T_tot * N_pre
  R = T_tot * R_pre 


  for counter in range(num*(skip+1)):
    Qnew =  np.zeros((m,m)) #np.matrix(np.copy(T))
    for i in range(m): # rows
      for j in range(m): # columns
        Qnew[i,j] = Q[i,j]

    for i in range(m): # rows
      for j in range(m): # columns
        if i == j: continue
        if Qnew[i,j] == 0.0: continue
        if Qnew[j,j] == 0.0: continue
        delta = random.expovariate(1.0/(Qnew[i,j])) - Qnew[i,j] # so that it's averaged at zero change, but has a minimum value of changing Q[j,j] down only to zero

        if np.isinf(delta): continue
        r = random.random()

        # NOTE: all this math is being done in logarithmic form first (log-likelihood)
        new_ij = N[i,j] * log(Qnew[i,j] + delta) - ((Qnew[i,j] + delta) * R[i])
        old_ij = N[i,j] * log(Qnew[i,j]) - ((Qnew[i,j]) * R[i])
        p_acc = (new_ij - old_ij) # + (new_jj - old_jj)
        #print log(r), p_acc
        if log(r) <= p_acc: # this can be directly compared to the log of the random variable
          #print "accepted"
          Qnew[i,i] = Qnew[i,i] - delta
          Qnew[i,j] = Qnew[i,j] + delta
        #else:
          #print "reject"


    #print Qnew 
    if counter % (skip+1) == 0: # then save this result for later analysis
      Q_mats.append(Qnew)
    Q = Qnew
  return Q_mats

def plot_conv(N_conv, R_conv, k_conv, conv_intervals, k_cell_conv, p_equil_conv):
  fig= plt.figure()
  cm = plt.get_cmap('tab20')
  ax = fig.add_subplot(5,1,1,)
  NUM_COLORS=20
  ax.set_color_cycle([cm(1.*j/NUM_COLORS) for j in range(NUM_COLORS)])
  for i in range(N_conv.shape[0]):
    for j in range(N_conv.shape[1]):
        if np.sum(N_conv[i][j][:]) != 0:
          label_string = 'Src: '+str(i) +',' + 'Dest: '+str(j)
          ax.plot(np.multiply(conv_intervals,2e-6), N_conv[i][j][:], label = label_string ,linestyle='-', marker="o", markersize = 1)
  ax.set_ylabel('N/T')
  ax.set_xlabel('time (ns)')
  #ax_settitle('Transition Count Convergence')
  box = ax.get_position()
  ax.set_position([box.x0,box.y0, box.width * 0.8, box.height])
  lgd1 = ax.legend(loc ='center left', bbox_to_anchor=(1, 0.5), ncol = 2)  
  ax.grid(b=True,axis = 'y', which = 'both')

  ax2 = fig.add_subplot(5,1,2, sharex = ax)
  ax2.set_color_cycle([cm(1.*j/NUM_COLORS) for j in range(NUM_COLORS)])
  for i in range(R_conv.shape[0]):
    for j in range(R_conv.shape[1]):
      if np.sum(R_conv[i][j][:]) != 0:
        label_string_2 = 'anchor ' +str(i) + ',' + 'Milestone '+str(j)
        ax2.plot(np.multiply(conv_intervals,2e-6), R_conv[i][j][:], label = label_string_2 ,linestyle='-', marker="o", markersize = 1)
  box = ax2.get_position()
  ax2.set_position([box.x0,box.y0, box.width * 0.8, box.height])
  lgd2 = ax2.legend(loc ='center left', bbox_to_anchor=(1, 0.5), ncol = 2)
  ax2.set_ylabel('R/T')
  ax2.set_xlabel('time (ns)')
  ax2.grid(b=True ,axis = 'y', which = 'both')
  #ax2_settitle('Incubation Time Convergence')

## k_cell and p_equil conv plots go here
 

  ax5 = fig.add_subplot(5,1,3, sharex = ax)
  ax5.set_color_cycle([cm(1.*j/NUM_COLORS) for j in range(NUM_COLORS)])
  for i in range(p_equil_conv.shape[0]):
        label_string_2 = 'anchor ' +str(i) 
        ax5.plot(np.multiply(conv_intervals,2e-6), p_equil_conv[i][:], label = label_string_2 ,linestyle='-', marker="o", markersize = 1)
  box = ax5.get_position()
  ax5.set_position([box.x0,box.y0, box.width * 0.8, box.height])
  lgd2 = ax5.legend(loc ='center left', bbox_to_anchor=(1, 0.5), ncol = 2)
  ax5.set_ylabel('p equil')
  ax5.set_xlabel('time (ns)')
  ax5.grid(b=True ,axis = 'y', which = 'both')

  ax4 = fig.add_subplot(5,1,4, sharex = ax)
  ax4.set_color_cycle([cm(1.*j/NUM_COLORS) for j in range(NUM_COLORS)])
  for i in range(k_cell_conv.shape[0]):
    for j in range(k_cell_conv.shape[1]):
      if np.sum(k_cell_conv[i][j][:]) != 0:
        label_string_2 = 'anchor ' +str(i) + ',' + 'Milestone '+str(j)
        ax4.plot(np.multiply(conv_intervals,2e-6), k_cell_conv[i][j][:], label = label_string_2 ,linestyle='-', marker="o", markersize = 1)
  box = ax4.get_position()
  ax4.set_position([box.x0,box.y0, box.width * 0.8, box.height])
  lgd2 = ax4.legend(loc ='center left', bbox_to_anchor=(1, 0.5), ncol = 2)
  ax4.set_ylabel('k cell')
  ax4.set_xlabel('time (ns)')
  ax4.grid(b=True ,axis = 'y', which = 'both')

  ax3 = fig.add_subplot(5,1,5, sharex = ax)
  ax3.plot(np.multiply(conv_intervals,2e-6), k_conv, linestyle='-', marker="o", markersize = 1)
  box = ax3.get_position()
  ax3.set_position([box.x0,box.y0, box.width * 0.8, box.height])
  #lgd3 = ax3.legend(loc ='center left', bbox_to_anchor=(1, 0.5), ncol = 2)
  ax3.set_ylabel('k off (s^-1)')
  ax3.set_xlabel('time (ns)')

  



  plt.savefig('convergence.png', bbox_extra_artists=(lgd1,lgd2), bbox_inches='tight', format='png', dpi=300)
  pickle.dump(fig, open('conv.fig.pickle', 'wb'))
  plt.show()
  return     


def main():
  parser = argparse.ArgumentParser(description='Analyzes SEEKR forward phase output to construct milestoning calculations on simulation analysis.')
  #parser.add_argument('-s', '--sinkrad', dest="sinkrad", nargs=1, help="Sink radius for the Markov Model", action=Once) # make sure this is only specified once, cause its easy to mix this one up with '--skip'
  onoff_group = parser.add_mutually_exclusive_group()
  onoff_group.add_argument('--on', dest="on", default=True, help="Perform a k-on calculation - where a low state is the sink", action="store_true")
  onoff_group.add_argument('--off', dest="off", default=False, help="Perform a k-off calculation - where a high state is the sink", action="store_true")
  onoff_group.add_argument('--free_energy', dest="free_energy", default=False, help="Calculate a free energy profile - where there are no sinks", action="store_true")
  parser.add_argument('--milestone_conv', dest='milestone_conv', default=False, help="If we are computing convergence properties for each milestone", action="store_true" )
  parser.add_argument('-m', '--milestones', dest="milestones", type=str, help="Milestones file") # This should contain most of what the user needs
  parser.add_argument('-b', '--bound_states', dest="bound_states", type=str, default="0", help="The milestone index of the bound state(s). If different bound states exist for different sites, then separate with a colon. For multiple bound states, separate with commas. Examples: '0', '1:2', '0:1,1:3'.")
  # TODO: escape state?
  parser.add_argument('--nobd', dest="nobd", help="Do not include the BD statistics in the calculation")
  parser.add_argument('--nomd', dest="nomd", help="Do not include the MD statistics in the calculation")
  parser.add_argument('--test', dest="test", help="Run the unittests", action="store_true")
  parser.add_argument('--doing_error', dest="doing_error", default=True, action = "store_true", help= "whether MCMC error estimate should be performed")
  parser.add_argument('--skip', dest='error_skip', type=int, default=1, help="define how many matrix samples to skip when doing MC error estimation" )
  parser.add_argument('--number', dest='error_number', type=int, default=1000, help="define how many matrices to sample for MC error estimation" )
  parser.add_argument('-v','--verbose', dest="verbose", help="verbose output", action="store_true")
  parser.add_argument('-i','--info', dest="info", help="Print information on transitions without computing anything, including the number of transition statistics read.", action="store_true")
  parser.add_argument('--conv_filename', dest='conv_filename', type=str, default="", help="If provided, an analysis of convergence of the computed value will be performed, and the results will be written to the specified file." )
  parser.add_argument('--conv_stride', dest='conv_stride', type=int, default=100, help="The stride through umbrella sampling statistics when a convergence analysis is being performed." )
  parser.add_argument('--rand_conv', dest='rand_conv', default=False,action="store_true", help="If we are computing randomized convergence properties" )
  parser.add_argument('--rand_conv_samples', dest='rand_conv_samples', type= int, default=0, help="number of random samples to take at each convergence interval" )
  parser.add_argument('--plt_name', dest='plt_name', type= str, default='', help="system name for use in plot titles" )


  args = parser.parse_args()
  args = vars(args) # convert to a dictionary
  #print "args:", args
  milestone_filename = args['milestones'] # parse the milestoning XML file
  doing_error = args['doing_error']
  verbose = False
  info_mode = False
  if args['verbose']:
    print "verbose mode activated."
    verbose = True

  if args['info']:
    print "'info' mode active, Only printing read information, not running any kinetics calculations."
    info_mode = True

  if args['test']: # if we are just going to run the unittests, then run and exit
    print "running Unittests..."
    runner = unittest.TextTestRunner()
    itersuite = unittest.TestLoader().loadTestsFromTestCase(Test_analyze) # because for some reason, unittest.main() doesn't work right
    runner.run(itersuite)
    exit()

  assert milestone_filename, "A milestones XML file must be provided"
  # figure out whether we are doing a k-on, k-off, or free energy profile calculation
  bd_time = 0.0 # something to allow the calculations to work
  error_skip = args['error_skip']
  error_number = args['error_number']
  conv_filename = args['conv_filename']
  conv_stride = args['conv_stride']
  milestone_conv= args['milestone_conv']
  #rand_conv= args['rand_conv']

  if args['off']:
    print "Running k-off calculations."
    calc_type = "off"

  elif args['free_energy']:
    print "Running free energy profile calculations."
    calc_type = "free_energy"
    assert not conv_filename, "Alert: free energy profile not currently supported in convergence mode. Please leave --conv_filename argument blank. Aborting..."

  #elif args['milestone_conv']:
  #  calc_type= "off"
  #rand_conv= args['rand_conv']

  else:
    print "Running k-on calculations."
    calc_type = "on" # by default
    bd_time = 1.0
  bound_dict = parse_bound_state_args(args['bound_states'])

  # open milestoning file and parse everything into a full milestoning model
  model = parse_milestoning_file(milestone_filename)


  #print len(model.sites[0].anchors)
  max_steps = read_transition_statistics_from_files(model, verbose)
  #max_steps = info_dict['max_steps']
  #print  len(model.sites[0].anchors[0].transitions)
  #if info_mode:
  #  print "Highest umbrella sampling step", info_dict['max_umbrella_steps']
  #  print "more info to be added later..."
  #  print "Aborting."
  #  exit()



  #print milestone_conv
  #print max_steps
  # starting kinetics analysis
  if calc_type == "off":
    p_equil, N, R, T, T_tot, Q, n_conv, r_conv, k_cell= analyze_kinetics(calc_type, model, bound_dict, doing_error=True, verbose=verbose, bd_time=bd_time, error_number=error_number, error_skip=error_skip)

    

    if doing_error == True:
      mfpt_list = []
      k_off_list = []
      #Q_mats = monte_carlo_milestoning_error(Q, N, R, p_equil,T_tot, num = error_number, skip = error_skip)
      Q_mats = monte_carlo_milestoning_nonreversible_error(Q, N, R, p_equil, T_tot, num = error_number, skip = error_skip)
      #print Q_mats
      for Q_err in Q_mats:
        T_err = calc_MFPT_vec(Q_err)
        #print "T", T_err
        mfpt_list.append(T_err[0])
        k_off_list.append(1e15/T_err[0])
      mfpt_std = np.std(mfpt_list)
      k_off_std = np.std(k_off_list)

    MFPT = T[0]
    k_off = 1e15/MFPT
    #print k_off_list
    print "avg k off", np.average(k_off_list)
    if doing_error == True:
      print "k_off: " , k_off," +- ", k_off_std, " s^-1" 
    else:
      print "k_off: " , k_off, " s^-1" 


  if milestone_conv == True:
    k_conv_file = open('k_off_conv.txt', 'w')
    conv_intervals = np.arange(conv_stride, max_steps, conv_stride)
    print conv_intervals
    print max_steps
    N_conv = np.zeros((15,15,len(conv_intervals)))
    R_conv = np.zeros((15,15,len(conv_intervals)))
    k_conv = np.zeros(len(conv_intervals))
    k_cell_conv = np.zeros((15,15,len(conv_intervals)))
    p_equil_conv = np.zeros((15,len(conv_intervals)))
    for interval_index in range(len(conv_intervals)):
      p_equil, N, R, T, T_tot, Q, n_conv, r_conv, k_cell = analyze_kinetics(calc_type, model, bound_dict, doing_error=False, verbose=verbose, max_steps=conv_intervals[interval_index], bd_time=bd_time, error_number=error_number, error_skip=error_skip)
     
      print n_conv 
      MFPT = T[0]
      k_off = 1e15/MFPT
      #print conv_intervals[interval_index], ",   ", k_off

      for index, x in np.ndenumerate(n_conv):
          N_conv[index[0]][index[1]][interval_index]=x
      for index2,y in np.ndenumerate(r_conv):
        R_conv[index2[0]][index2[1]][interval_index]= y 
      for index3,z in np.ndenumerate(k_cell):
        k_cell_conv[index3[0]][index3[1]][interval_index]= z
      for index4,j in np.ndenumerate(p_equil):
        p_equil_conv[index4[0]][interval_index]= j   
      k_conv[interval_index]=k_off    
      

      if verbose: print p_equil_conv.shape
      
      #k_conv_file.write(str(interval)+'\t'+(str(k_off))+'\n')

      #print np.shape(n_conv)
      #pprint(n_conv)
      #print np.shape(r_conv)
      #np.concatenate((N_conv, n_conv), )
      #np.concatenate((R_conv, r_conv), )
      #k_conv.append(k_off)
    #k_conv_file.close()
    #pprint(N_conv)
    #pprint(R_conv)
    MFPT = T[0]
    k_off = 1e15/MFPT
    print "k_off: " , k_off, "s^-1"

    print 'conv intervals' , np.multiply(conv_intervals,2e-6)
    #print N_conv
    plot_conv(N_conv, R_conv, k_conv, conv_intervals, k_cell_conv, p_equil_conv)
  #else:
  #  p_equil, N, R, T, n_conv, r_conv, k_cell, p_equil= analyze_kinetics(calc_type, model, bound_dict, doing_error=False, verbose=verbose, bd_time=bd_time, error_number=error_number, error_skip=error_skip)

   #  MFPT = T[0]
   #  k_off = 1e15/MFPT
   #  print "k_off: " , k_off, "s^-1"
   # # print max_steps, ",   ", k_off


if __name__ == "__main__": main()

