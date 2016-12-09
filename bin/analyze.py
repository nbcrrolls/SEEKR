#!/usr/bin/python
# encoding: utf-8
'''
analyze.py

analyzes forward stage output to construct markov model for simulation analysis

analyze -- combines the statistics generated from MD and BD milestoning to estimate kinetic rates

@author:     Lane Votapka

@copyright:  2015 Amarolab UCSD. All rights reserved.

@license:    Licensed under the Apache License 2.0 - http://www.apache.org/licenses/LICENSE-2.0

@contact:    lvotapka100@gmail.com
@deffield    updated: Feb 6, 2015
'''

import os, sys, re, random, unittest, glob, tempfile
import argparse
from string import Template
from math import exp, log
from pprint import pprint
import numpy as np
from scipy import linalg
from copy import deepcopy
import xml.etree.ElementTree as ET
from subprocess import check_output

k_boltz = 1.3806488e-23
R_GAS_CONSTANT = 0.0019872 # in kcal/mol*K
DEFAULT_TEMP = 0.0 # temperature to assign if not in the xml file
DEFAULT_MD_TIME_FACTOR = 2.0
DEFAULT_BD_TIME_FACTOR = 1000.0
DEFAULT_SHAPE = "sphere"
#GRID_TRANSITION = re.compile("GRID TRANSITION")
GRID_TRANSITION = "SEEKR: GRID TRANSITION: "
GRID_TRANSITION_COMPILE = re.compile(GRID_TRANSITION)

FORWARD_OUTPUT_GLOB = "fwd_rev*.out.*" # used to find all the forward phase simulation output

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
    self.num_milestones += site.num_milestones
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

  #def construct_transition_matrices
  #...



class Site():
  ''' object representing a site: a collection of related milestones that share, say, the same origin '''
  def __init__(self, index, name):
    self.milestones = []
    self.num_milestones = 0
    self.index = index
    self.name = name
    #self.center ?
    #..

  def add_milestone(self, milestone):
    self.milestones.append(milestone)
    self.num_milestones += 1

class Milestone():
  ''' object representing a single milestone '''
  def __init__(self, index, shape, end, md, bd, fullname, directory, siteindex, sitename, anchor=None, normal=None, radius=None, ):
    self.index = index # most of this information is read from the provided milestone.xml file
    self.shape = shape
    self.end = end.lower()
    self.fullname = fullname
    self.directory = directory
    self.anchor = anchor
    self.normal = normal
    self.radius = radius
    self.bd = bd
    self.md = md
    self.site = siteindex
    self.sitename = sitename
    self.transitions = [] # all the transition statistics associated with this milestone


  def parse_md_transitions(self):
    'find all forward phase simulation output files'
    forward_dir_glob = os.path.join(self.directory,'md','fwd_rev',FORWARD_OUTPUT_GLOB)
    #print forward_dir_glob
    forward_output_filenames = glob.glob(forward_dir_glob)
    # read files and sift out the transition lines
    transition_lines = [] # a list containing all lines of the output that mark transitions
    for filename in forward_output_filenames:
     # print 'output filename', filename
      for line in open(filename,'r'):
        if re.match(GRID_TRANSITION_COMPILE, line):
          transition_lines.append(line)
      #pprint(transition_lines)
    # feed the lines into the Transition object
    for line in transition_lines:
      self.transitions.append(Transition(line))
      
    return

  def get_md_transition_statistics(self, md_time_factor=DEFAULT_MD_TIME_FACTOR):
    'parse the transition data to obtain transition counts'
    counts = {} # all the sources and their destinations
    total_counts = {} # keeps track of all counts to any destination
    total_times = {} # the total time of all counts, to be averaged later
    avg_times = {} # the times to transition out of each source
    site_index = self.site
    site_name = self.sitename
    for transition in self.transitions:
      source = transition.src
      dest = transition.dest
      time = transition.time
      #src_key = '%d_%d' % (site_index, source)
      src_key = '%s_%d' % (site_name, source)
      #dest_key = '%d_%d' % (site_index, dest)
      dest_key = '%s_%d' % (site_name, dest)
      if src_key in counts.keys():
        if dest_key in counts[src_key].keys():
          counts[src_key][dest_key] += 1
        else:
          counts[src_key][dest_key] = 1
        total_times[src_key] += time * md_time_factor
        total_counts[src_key] += 1
      else:
        counts[src_key] = {dest_key:1}
        total_counts[src_key] = 1
        total_times[src_key] = time * md_time_factor

    for src_key in total_times.keys():
      avg_times[src_key] = total_times[src_key] / total_counts[src_key]

    return counts, total_counts, total_times, avg_times
    
  def get_bd_transition_statistics(self, results_filename=os.path.join("bd","results.xml"), bd_time=0.0):
    'read the BD results.xml file for the anchors to determine transition statistics and times for the BD stage'
    bd_results_filename = os.path.join(self.directory, results_filename)
    counts = {}
    total_counts = {}
    total_times = {}
    avg_times = {}
    src_key = '%s_%s' % (self.sitename, self.index)
    counts[src_key], total_times[src_key] = parse_bd_results(bd_results_filename)
    total_counts[src_key] = 0
    for dest_key in counts[src_key].keys():
      total_counts[src_key] += counts[src_key][dest_key]
    if total_times[src_key] == None:
      avg_times[src_key] = bd_time
      total_times[src_key] = bd_time
    else:
      avg_times[src_key] = total_times[src_key] / total_counts[src_key]
    
    return counts, total_counts, total_times, avg_times

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
    self.ligand_com = linedict['ligand COM'].strip().split()
    self.receptor_com = linedict['receptor COM'].strip().split()
    self.receptor_start_com = linedict['receptor start COM'].strip().split()

  def print_status(self):
    print "src:", self.src
    print "dest:", self.dest
    print "cur_step:", self.cur_step
    print "time:", self.time
    print "ligand_com:", self.ligand_com
    print "receptor_com:", self.receptor_com
    print "receptor_start_com:", self.receptor_start_com

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
    for milestone in branch: # iterate thru all milestones
      if milestone.tag != "milestone": continue
      index = milestone.find('index').text.strip()
      shape_xml = milestone.find('shape')
      if shape_xml != None:
        shape = shape_xml.text.strip()
      else:
        shape = DEFAULT_SHAPE # by default
      anchor = milestone.find('anchor').text.strip()
      fullname = milestone.find('fullname') # this may or may not exist in old milestoning files
      if fullname != None:
        fullname = fullname.text.strip() # parse all the attributes of the milestone
      else: # then just name it by its anchor
        fullname = str(anchor)
      directory_text = milestone.find('directory').text
      if directory_text:
        directory = directory_text.strip() # directory where this milestone is located
      else:
        directory = None
      radius = None
      normal = None
      end = milestone.find('end').text.strip()
      bd =boolean(milestone.find('bd').text.strip())
      #if milestone.find('md'):
      md =boolean( milestone.find('md').text.strip())
      #print 'MD?' , md
      #else:
        #md = "False"
      if shape == "plane": # get other information based on the milestone shape
        anchor = milestone.find('anchor').text.strip()
        normal = milestone.find('normal').text.strip()
      elif shape == "sphere":
        anchor = milestone.find('anchor').text.strip()
        radius = milestone.find('radius').text.strip()
      elif shape == "rotational":
        pass
      milestone_obj = Milestone(index, shape, end, md, bd, fullname, directory, site_counter, site_name, anchor, normal, radius) # create the Milestone object
      site_obj.add_milestone(milestone_obj) # add milestone to the site
    model.add_site(site_obj) # add site to the model
    site_counter += 1
    
    # parse the b-surface
    model.b_surface_milestone = Milestone(index="0", shape="sphere", end="True", md=False, bd=True, fullname="b_surface", directory="b_surface", siteindex=0, sitename="b_surface") # create the Milestone object
    
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


def trans_dict_to_matrix(trans_dict):
  '''given a transition (or count) dictionary, will return a transition (or count) matrix'''
  index_dict = {}
  trans_dict_keys = trans_dict.keys()
  trans_dict_keys = [key.replace('inf', 'inf_0') for key in trans_dict_keys]
  n = len(trans_dict_keys)
  i = 0
  trans_dict_keys = sorted(trans_dict_keys, key=lambda keystring:keystring.split('_')[0]) # first sort by site
  trans_dict_keys = sorted(trans_dict_keys, key=lambda keystring:int(keystring.split('_')[1])) # then sort by milestone
  #count_dict_keys = sorted(count_dict.keys(), key=lambda keystring:keystring )#int(keystring.split('_')[1])) # then sort by milestone
  #count_dict_keys.sort()
  for key in trans_dict_keys:
    if key == 'inf_0': key = 'inf'
    index_dict[i] = key
    i += 1

  #print "index dict:", index_dict

  trans_matrix = np.zeros((n,n))


  for i in range(n):
    for f in range(n):
      if index_dict[f] in trans_dict[index_dict[i]].keys():
        trans = trans_dict[index_dict[i]][index_dict[f]]
        trans_matrix[f,i] = trans
        #sum_vector[f] += count
      else:
        trans_matrix[f,i] = 0.0

  '''for i in range(n):
    for f in range(n):
      trans_matrix[i][f] = count_matrix[i][f] / sum_vector[f]'''

  return trans_matrix, index_dict
  
def trans_dict_to_q0_vector(trans_dict, index_dict):
  '''given a transition matrix, will return a vector of starting fluxes based on b-surface stats.'''
  n = len(index_dict.keys())
  src_key = trans_dict.keys()[0] # the key that refers to the b-surface trans dict
  q0_vector = np.zeros((n,1))
  for i in range(n): # move down the vector element by element
    dest_key = index_dict[i]
    if dest_key in trans_dict[src_key].keys():
      q0_vector[i,0] = trans_dict[src_key][dest_key]
  return q0_vector
  

def avg_t_vector(time_dict, index_dict):
  ''' given a dictionary, will construct a vector in the correct order'''
  n = len(index_dict.keys())
  t_count = []
  for i in range(n):
    if index_dict[i] in time_dict.keys():
      t_count.append(time_dict[index_dict[i]])
    else:
      t_count.append(0)
  t_matrix = np.array([t_count]).T #* (1/float(sum_total))
  #print "p_matrix", p_matrix
  return np.matrix(t_matrix)
  
def parse_bound_state_args(bound_args):
  bound_dict = {}
  bound_pairs = bound_args.split(',')
  for pair in bound_pairs:
    site_index = pair.split(':')
    if len(site_index) == 1:
      site = 'all'
      index = site_index[0]
    elif len(site_index) == 2:
      site = site_index[0]
      index = site_index[1]
    bound_dict[site] = index

  return bound_dict

def run_compute_rate_constant(results_filename, browndye_bin_dir=""):
  'runs the Browndye program compute_rate_constant to find the value k(b)'
  process_trajectories = os.path.join(browndye_bin_dir, "compute_rate_constant")
  cmd = "%s < %s" % (process_trajectories, results_filename)
  output_string = check_output(cmd, shell=True) # run the Browndye process command
  root = ET.fromstring(output_string) # read the XML string
  rate = root[0] # the first <rate> tag, because it doesn't really matter which one we choose
  rate_constant = rate.find('rate-constant')
  k_on_tag = rate_constant.find('mean')
  reaction_probability = rate.find('reaction-probability')
  beta_tag = reaction_probability.find('mean')
  assert (k_on_tag != None) and (beta_tag != None), "Alert: Unable to find rate constant <mean> tags after running compute_rate_constant on file %s" % results_filename
  k = float(k_on_tag.text)
  beta = float(beta_tag.text)
  k_b = k / beta # the flux to the b-surface
  return k_b

def count_mat_to_rate_mat(N, avg_t):
  ''' converts a count matrix N to a markov rate matrix Q (where each entry is an effective kinetic rate constant)'''
  n = np.shape(N)[0]
  Q = np.matrix(np.zeros((n,n)))
  sum_vector = np.zeros(n)
  for i in range(n): # first make the sum of all the counts
    for j in range(n):
      sum_vector[j] += N[i,j]


  for i in range(n):
    for j in range(n):
      if j == i: continue
      if sum_vector[j] == 0 or avg_t[j] == 0.0:
        Q[i,j] = 0.0
      else:
        Q[i,j] = N[i,j] / (sum_vector[j] * avg_t[j])
      Q[j,j] -= Q[i,j]

  return Q, sum_vector

def rate_mat_to_prob_mat(Q, calc_type='on'):
  ''' converts a rate matrix Q into probability matrix (kernel) K and an incubation time vector'''
  n = Q.shape[0]
  P = np.matrix(np.zeros((n,n)))
  K = np.matrix(np.zeros((n,n)))
  sum_vector = np.zeros(n)
  avg_t = np.zeros(n)
  for i in range(n): # first make the sum of all the rates
    for j in range(n):
      if j == i: continue
      sum_vector[j] += Q[i,j]


  for i in range(n):
    for j in range(n):
      if j == i: continue
      if sum_vector[j] == 0:
        K[i,j] = 0.0
        #avg_t[j] = 0.0
      else:
        K[i,j] = Q[i,j] / sum_vector[j]

    if sum_vector[i] != 0.0:
      avg_t[i] = 1.0 / sum_vector[i]
    else: # then the sum vector goes to zero, make the state go to itself
      if calc_type == "on":
        K[i,i] = 1.0
      elif calc_type == "off":
        K[i,i] = 0.0
      #avg_t[i] = something???

  return K, avg_t

def monte_carlo_milestoning_nonreversible_error(N, avg_t, num = 20, skip = 0):
  ''' Samples a distribution of rate matrices that are nonreversible
      using Markov Chain Monte Carlo method

      The distribution being sampled is:
      p(Q|N) = p(Q)p(N|Q)/p(N) = p(Q) PI(q_ij**N_ij * exp(-q_ij * N_i * t_i))

      N = count matrix
      avg_t = incubation times

  '''
  Q0, N_sum = count_mat_to_rate_mat(N, avg_t) # get a rate matrix and a sum of counts vector
  m = N.shape[0] # the size of the matrix
  Q = Q0
  Q_mats = []

  for counter in range(num*(skip+1)):
    Qnew =  np.matrix(np.zeros((m,m))) #np.matrix(np.copy(T))
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
        new_ij = N[i,j] * log(Qnew[i,j] + delta) - ((Qnew[i,j] + delta) * N_sum[j] * avg_t[j])
        old_ij = N[i,j] * log(Qnew[i,j]) - ((Qnew[i,j]) * N_sum[j] * avg_t[j])
        p_acc = (new_ij - old_ij) # + (new_jj - old_jj)
        if log(r) <= p_acc: # this can be directly compared to the log of the random variable
          Qnew[j,j] = Qnew[j,j] - delta
          Qnew[i,j] = Qnew[i,j] + delta

    if counter % (skip+1) == 0: # then save this result for later analysis
      Q_mats.append(Qnew)
    Q = Qnew
  return Q_mats

def get_beta_from_K_q0(K, q0, bound_indices):
  'given a transition matrix K and starting vector q_0, returns a beta value for k-ons'
  K_inf = np.matrix(K) ** 99999999
  #n,m = K_inf.shape
  q_inf = np.dot(K_inf, q0)
  #print "q_inf:", q_inf
  beta = 0.0
  for bound_index in bound_indices:
    beta += q_inf[bound_index, 0]
  return beta

def make_free_energy_profile_boundaries(K, bound_indices, inf_index):
  n = K.shape[0]
  for bound_index in bound_indices:
    K[bound_index,bound_index] = 0.0001
  
  for i in range(n):
    col_sum = 0.0
    
    for j in range(n):
      if j == inf_index:
        K[j, i] = 0.0
      col_sum += K[j, i]
    
    if col_sum > 0.0:
      K[:,i] /= col_sum
    
  return K
  

def main():
  parser = argparse.ArgumentParser(description='Analyzes SEEKR forward phase output to construct milestoning calculations on simulation analysis.')
  #parser.add_argument('-s', '--sinkrad', dest="sinkrad", nargs=1, help="Sink radius for the Markov Model", action=Once) # make sure this is only specified once, cause its easy to mix this one up with '--skip'
  onoff_group = parser.add_mutually_exclusive_group()
  onoff_group.add_argument('--on', dest="on", default=True, help="Perform a k-on calculation - where a low state is the sink", action="store_true")
  onoff_group.add_argument('--off', dest="off", default=False, help="Perform a k-off calculation - where a high state is the sink", action="store_true")
  onoff_group.add_argument('--free_energy', dest="free_energy", default=False, help="Calculate a free energy profile - where there are no sinks", action="store_true")
  parser.add_argument('-m', '--milestones', dest="milestones", type=str, help="Milestones file") # This should contain most of what the user needs
  parser.add_argument('-b', '--bound_states', dest="bound_states", type=str, default="0", help="The milestone index of the bound state(s). If different bound states exist for different sites, then separate with a colon. For multiple bound states, separate with commas. Examples: '0', '1:2', '0:1,1:3'.")
  # TODO: escape state?
  parser.add_argument('--nobd', dest="nobd", help="Do not include the BD statistics in the calculation")
  parser.add_argument('--nomd', dest="nomd", help="Do not include the MD statistics in the calculation")
  parser.add_argument('--test', dest="test", help="Run the unittests", action="store_true")
  parser.add_argument('--skip', dest='error_skip', type=int, default=1, help="define how many matrix samples to skip when doing MC error estimation" )
  parser.add_argument('--number', dest='error_number', type=int, default=1000, help="define how many matrices to sample for MC error estimation" )
  parser.add_argument('-v','--verbose', dest="verbose", help="verbose output", action="store_true")

  args = parser.parse_args()
  args = vars(args) # convert to a dictionary
  print "args:", args
  milestone_filename = args['milestones'] # parse the milestoning XML file
  verbose = False
  if args['verbose']:
    print "verbose mode activated."
    verbose = True
    
  if args['test']: # if we are just going to run the unittests, then run and exit
    print "running Unittests..."
    runner = unittest.TextTestRunner()
    itersuite = unittest.TestLoader().loadTestsFromTestCase(Test_md) # because for some reason, unittest.main() doesn't work right
    runner.run(itersuite)
  assert milestone_filename, "A milestones XML file must be provided"
  # figure out whether we are doing a k-on, k-off, or free energy profile calculation
  bd_time = 0.0 # something to allow the calculations to work
  if args['off']:
    print "Running k-off calculations."
    calc_type = "off"
    
  elif args['free_energy']:
    print "Running free energy profile calculations."
    calc_type = "free_energy"
  else:
    print "Running k-on calculations."
    calc_type = "on" # by default
    bd_time = 1.0
    
  error_skip = args['error_skip']
  error_number = args['error_number']
    
  bound_dict = parse_bound_state_args(args['bound_states'])

  # open milestoning file and parse everything into a full milestoning model
  model = parse_milestoning_file(milestone_filename)
  #model.make_directories()
  # gather MD data - fill into model

  counts = {}; total_counts = {}; total_times = {}; avg_times = {}; trans = {}
  end_indeces = []
  for site in model.sites:
    for milestone in site.milestones:
     # print "directory:", milestone.directory
     # print 'MD2?', milestone.md
     # print 'BD2?', milestone.bd
      if milestone.md == True and milestone.directory:
        print 'parsing md transitions for:Anchor', milestone.fullname
        milestone.parse_md_transitions()
        this_counts, this_total_counts, this_total_times, this_avg_times = milestone.get_md_transition_statistics(model.md_time_factor) # find the count statistics
        total_counts = add_dictionaries(total_counts, this_total_counts)
        total_times = add_dictionaries(total_times, this_total_times)
        for src_key in this_counts.keys():
          if src_key in counts.keys():
            counts[src_key] = add_dictionaries(counts[src_key], this_counts[src_key])
          else:
            counts[src_key] = this_counts[src_key]
        #print "len(transitions)", len(milestone.transitions)
      if milestone.end == "true": # then its an end milestone, and there will be no transitions out of it
        end_indeces.append('%d_%d' % (milestone.site, int(milestone.index)))
        
      if milestone.bd == True and milestone.directory:
        this_counts, this_total_counts, this_total_times, this_avg_times = milestone.get_bd_transition_statistics(bd_time=bd_time)
        total_counts = add_dictionaries(total_counts, this_total_counts)
        total_times = add_dictionaries(total_times, this_total_times)
        for src_key in this_counts.keys():
          if src_key in counts.keys():
            counts[src_key] = add_dictionaries(counts[src_key], this_counts[src_key])
          else:
            counts[src_key] = this_counts[src_key]
        
  for src_key in total_times.keys(): # construct the average incubation times
    avg_times[src_key] = total_times[src_key] / total_counts[src_key]

  for src_key in counts.keys():
    temp = {}
    for dest_key in counts[src_key].keys(): # make the transition probability dictionary
      temp[dest_key] = float(counts[src_key][dest_key]) / float(total_counts[src_key])
      if dest_key in end_indeces: # if this milestone is going toward an end index, then create it in the count matrix and set the time to zero, transitioning to the one it came from
        #trans[dest_key] = {src_key: 0.9, dest_key: 0.1}
        #avg_times[dest_key] = 0.0
        pass
    trans[src_key] = temp
  #print "end_indeces:", end_indeces
  
  
  
  # b-surface milestone
  b_surface_counts, b_surface_total_counts, b_surface_total_times, b_surface_avg_times = model.b_surface_milestone.get_bd_transition_statistics("results.xml", bd_time=bd_time)
  src_key = b_surface_counts.keys()[0]
  b_surface_trans = {src_key:{}}
  bound_indices = []
  bound_keys = []
  for dest_key in b_surface_counts[src_key].keys():
    b_surface_trans[src_key][dest_key] = float(b_surface_counts[src_key][dest_key]) / float(b_surface_total_counts[src_key])
  
  
  radius_dict = {}
  
  if calc_type == "on": # TODO: these need to be added to the counts matrix somehow for the MC error
    trans['inf'] = {'inf':1.0}
    counts['inf'] = {'inf':1e99}
    avg_times['inf'] = 1.0
    found_sink = False
    for site in model.sites:
      for milestone in site.milestones:
        key = "%s_%s" % (site.name, milestone.index)
        if site.name in bound_dict.keys():
          if milestone.index in bound_dict[site.name]: # then make it the bound state by altering the transition properties.
            bound_keys.append(key)
            trans[key] = {key:1.0}
            counts[key] = {key:1e99}
            avg_times[key] = 1.0
            found_sink = True
            if verbose: print "setting site: %s to be the sink state" % key
        if 'all' in bound_dict.keys():
          if milestone.index in bound_dict['all']:
            #key = "%s_%s" % (site.name, milestone.index)
            bound_keys.append(key)
            trans[key] = {key:1.0}
            counts[key] = {key:1e99}
            avg_times[key] = 1.0
            found_sink = True
            if verbose: print "setting site: %s to be the sink state" % key
        
        if milestone.shape == "sphere":
          radius_dict[key] = milestone.radius
      
    assert found_sink, "Alert: no bound state found matching criteria: %s. Make sure that the site names and milestone indices match." % args['bound_states']
  elif calc_type == "off":
    trans['inf'] = {'inf':0.0}
    counts['inf'] = {'inf':0}
    # TODO: fill out more here
    avg_times['inf'] = 0.0
    
  elif calc_type == "free_energy": # then the user wants a free energy profile
    trans['inf'] = {'inf':0.0}
    counts['inf'] = {'inf':0}
    avg_times['inf'] = 1.0
    found_sink = False
    for site in model.sites:
      for milestone in site.milestones:
        key = "%s_%s" % (site.name, milestone.index)
        if site.name in bound_dict.keys():
          
          if 'inf' in trans[key].keys(): # we don't want anything transitioning to the infinite state
            trans[key]['inf'] = 0.0
          bound_keys.append(key)
          
        if 'all' in bound_dict.keys():
          if milestone.index in bound_dict['all']:
            #key = "%s_%s" % (site.name, milestone.index)
            bound_keys.append(key)
            
        if milestone.shape == "sphere":
          radius_dict[key] = milestone.radius
    
    '''
    found_sink = False
    for site in model.sites:
      for milestone in site.milestones:
        if site.name in bound_dict.keys():
          if milestone.index in bound_dict[site.name]: # then make it the bound state by altering the transition properties.
            key = "%s_%s" % (site.name, milestone.index)
            bound_keys.append(key)
            #trans[key] = {key:0.0} # infinity is now a sink state
            #counts[key] = {key:0}
            #avg_times[key] = 0.0
            found_sink = True
            print "setting site: %s to be the sink state" % key
        if 'all' in bound_dict.keys():
          if milestone.index in bound_dict['all']:
            key = "%s_%s" % (site.name, milestone.index)
            bound_keys.append(key)
            #trans[key] = {key:0.0}
            #counts[key] = {key:0}
            #avg_times[key] = 0.0
            found_sink = True
            print "setting site: %s to be the sink state" % key
    assert found_sink, "Alert: no bound state found matching criteria: %s. Make sure that the site names and milestone indices match." % args['bound_states']
  '''
  if verbose:
    print "bound_dict:", bound_dict
    print "trans:"
    print trans
  
    print "radius_dict:"
    print radius_dict
  
  #avg_t = np.zeros((n,1))
  K, index_dict = trans_dict_to_matrix(trans)
  n,m = K.shape
  N, N_index_dict = trans_dict_to_matrix(counts)
  avg_t = avg_t_vector(avg_times, index_dict)
  inf_index = None
  for key in index_dict.keys():
    if index_dict[key] == 'inf':
      assert inf_index == None, "cannot have more than one infinite state."
      inf_index = key
  
  if verbose:
    print "n:", n, "m:", m
    print "K:", K
    print "index_dict:", index_dict
    print "N:", 
    print N
    print "N_index_dict:", N_index_dict
    print "counts:", counts
    print "inf_index:", inf_index
  
  # Now sampling MC matrices
  print "Now sampling rate matrices for MC error estimation."
  Q_mats = monte_carlo_milestoning_nonreversible_error(N, avg_t, num = error_number, skip = error_skip)
  
  
  for index_key in index_dict.keys():
    if index_dict[index_key] in bound_keys:
      bound_indices.append(index_key)
  
  if calc_type == "on": # then get the starting probabilities from the b-surface statistics
    q0 = trans_dict_to_q0_vector(b_surface_trans, index_dict)
  
  elif calc_type == "off" or calc_type == "free_energy":
    site_count = 0.0
    q0 = np.zeros((len(index_dict.keys()),1))
    for site in model.sites:
      for milestone in site.milestones:
        if site.name in bound_dict.keys():
          if milestone.index in bound_dict[site.name]: # then make it the bound state by altering the transition properties.
            key = "%s_%s" % (site.name, milestone.index)
            for i in index_dict.keys():
              if index_dict[i] == key:
                break
            q0[i,0] = 1.0
            site_count += 1.0
            if verbose: print "setting site: %s to be a starting state" % key
            
        if 'all' in bound_dict.keys():
          if milestone.index in bound_dict['all']:
            key = "%s_%s" % (site.name, milestone.index)
            for i in index_dict.keys():
              if index_dict[i] == key:
                break
            q0[i,0] = 1.0
            site_count += 1.0
            if verbose: print "setting site: %s to be a starting state" % key
    q0 = q0 / site_count # normalize
  
    
  
  # === NEED TO CLEAN UP BEYOND THIS POINT ===
  if verbose:
    print 'counts:' , counts
    print "avg_times dictionary:", avg_times
    print "avg_t vector", avg_t
    print "transitions:"
    pprint(K)
    print "index_dict:"
    pprint(index_dict)
    print "q0"
    print q0
    print "bound_keys:", bound_keys
    print "bound_indices:", bound_indices

  if calc_type == "on":
    beta = get_beta_from_K_q0(K, q0, bound_indices)
    betas = []
    for Q in Q_mats:
      #print "Q:", Q
      new_K, new_avg_t = rate_mat_to_prob_mat(Q)
      new_beta = get_beta_from_K_q0(new_K, q0, bound_indices)
      betas.append(new_beta)
      
    beta_std = np.std(betas)
    #print "len(betas):", len(betas)
    print "np.mean(betas):", np.mean(betas), "np.std(betas)", np.std(betas)

    print "beta:", beta, "+/-", beta_std
    k_b = run_compute_rate_constant(results_filename=os.path.join("b_surface", "results.xml"), browndye_bin_dir="")
    if verbose: print "k(b):", k_b
    k_on = k_b * beta
    k_on_std = k_b * beta_std
    print "k-on:", k_on, "+/-", k_on_std, "M^-1 s^-1"
    
    # now calculate p stationary
    #p = np.zeros((n,1))
    #print "p:"
    #total_p = 0.0
    #for i in range(n):
      #print q_inf[i,0]
      #p[i,0] = q_inf[i,0] * avg_t[i,0]
      #total_p += p[i,0]
  
    #print p
  '''
  for i in range(n):
    p[i,0] /= total_p
    print p[i,0]
  '''

  # now calculate MFPT
  #sink = end_indeces[-1] # just choose the last one to be a sink state
  #trans[sink] = {} # not going anywhere
  #t_mat_sink, index_dict = trans_dict_to_matrix(trans)
  #print "t_mat_sink:"
  #print t_mat_sink
  if calc_type == "off":
    #t_mat_sink = np.matrix(t_mat_sink)
    I = np.matrix(np.identity(n))
    aux = np.linalg.solve(I - K.T, avg_t)
    if verbose: print "aux:", aux
    mfpt = q0.T.dot(aux)
    koffs = []
    mfpts = []
    for Q in Q_mats:
      new_K, new_avg_t = rate_mat_to_prob_mat(Q, calc_type)
      new_mfpt = q0.T.dot(np.linalg.solve(I - new_K.T, new_avg_t))
      mfpts.append(new_mfpt)
      koffs.append(1e15/new_mfpt)
    
    #print "np.mean(koffs):", np.mean(koffs)
    print "MFPT:", mfpt, "+/-", np.std(mfpts), "fs"
    print "k-off:", 1e15/mfpt, "+/-", np.std(koffs), "s^-1"

  '''
  # now calculate flux quantity
  for sink in end_indeces: # all end states are sinks
    trans[sink] = {sink:1.0} # all sink states go to themselves
  t_mat_flux, index_dict = trans_dict_to_matrix(trans)

  t_mat_flux_inf = np.matrix(t_mat_flux) ** 2000000
  q_flux = t_mat_flux_inf* q0

  print "q_flux:"
  pprint(q_flux)
  print "q_flux[0]:", float(q_flux[0])
  print "q_flux[-1]:", float(q_flux[-1])

  print "trans:"
  pprint(trans)
  '''
  
  if calc_type == 'free_energy':
    pstat = np.matrix(np.zeros((n, 1)))
    delta_G = np.matrix(np.zeros((n,1)))
    
    K = make_free_energy_profile_boundaries(K, bound_indices, inf_index)
    K_inf = np.matrix(K) ** 99999999
    qstat = np.dot(K_inf,q0)
    
    #print "model.temperature:", model.temperature
    #print "type(model.temperature):", type(model.temperature)
    for i in range(n):
      pstat[i,0] = qstat[i,0] * avg_t[i,0]
      
    
    # TODO: change this hardcode to be a more general bound-state index instead of "1"
    pstat_ref = pstat[1,0]
    for i in range(n):
      if pstat[i,0] == 0.0:
        delta_G[i,0] = 1e99
      else:
        delta_G[i,0] = -model.temperature * R_GAS_CONSTANT * log(pstat[i,0] / pstat_ref) # in kcal/mol
    
    
    pstats = []
    delta_Gs = []
    for Q in Q_mats:
      new_K, new_avg_t = rate_mat_to_prob_mat(Q)
      
      new_K = make_free_energy_profile_boundaries(new_K, bound_indices, inf_index)
      new_pstat = np.matrix(np.zeros((n, 1)))
      new_delta_G = np.matrix(np.zeros((n, 1)))
      new_K_inf = np.matrix(new_K) ** 99999999
      new_qstat = np.dot(new_K_inf,q0)
      
      for i in range(n):
        new_pstat[i,0] = new_qstat[i,0] * new_avg_t[i]
      pstats.append(new_pstat)
      
      for i in range(n):
        if new_pstat[i,0] == 0.0:
          new_delta_G[i,0] = 1e99
        else:
          new_delta_G[i,0] = -model.temperature * R_GAS_CONSTANT * log(new_pstat[i,0] / pstat_ref) # in kcal/mol
      delta_Gs.append(new_delta_G)
         
    pstat_std = np.matrix(np.zeros((n, 1)))
    delta_G_std = np.matrix(np.zeros((n, 1)))
    for i in range(n):
      pstat_std_i = []
      delta_G_std_i = []
      for pstat in pstats:
        pstat_std_i.append(pstat[i,0])
      for mc_delta_G in delta_Gs:
        delta_G_std_i.append(mc_delta_G[i,0])
        
      pstat_std[i,0] = np.std(pstat_std_i)
      delta_G_std[i,0] = np.std(delta_G_std_i)
    
    if verbose:
      print "K:"
      pprint(K)
      print "K_inf:"
      pprint(K_inf)
      print "qstat:", qstat
      print "pstat:", pstat
      print "delta_G:", delta_G
      print "pstat_std:", pstat_std
      print "delta_G_std:", delta_G_std
    
    print "Free energy profile (kcals/mol):"
    print "radius:\tdelta G\t+/-"
    for i in range(n):
      radius = "???"
      if index_dict[i] in radius_dict.keys():
        radius = radius_dict[index_dict[i]]
      if delta_G[i,0] >= 1e99: continue # don't print huge numbers
      print '%s\t%2.3f\t%2.3f' % (radius, delta_G[i,0], delta_G_std[i,0])

  # ideas for analyze.py:
  # - lots of information in the milestones.xml file: directories, temperature, (check)
  # - option to choose bound state (check)
  # - lots of defaults and automation (check)
  # - classes for model, sites, milestones, etc.: all into a single object so that all information can be centralized before processing (check)
  # - unittests
  # - Warnings and alerts for unfinished forward trajectories



class Test_analyze(unittest.TestCase):
  # several test cases to ensure the functions in this module are working properly
  def test_main(self):
    print "Running unit tests for analyze.py"
    pass

  def test_parse_milestoning_file(self):
    test_milestoning_text = '''<?xml version="1.0" ?>
    <root><temperature>300</temperature><site><milestone><index>0</index> <fullname> milestone0 </fullname> <shape> plane </shape> <anchor> (0.0,0.0,-40.68) </anchor> <end> True </end> <bd> False </bd> <normal> (0.0,0.0,1.0) </normal> </milestone>
    <milestone><index>1</index> <fullname> milestone1 </fullname> <shape> sphere </shape> <anchor> (0.0,0.0,-30.68) </anchor> <end> False </end> <bd> False </bd> <radius> 3.1415 </radius> </milestone></site></root>
    '''
    test_milestoning_filename = '/usr/tmp/test_milestoning.xml'
    test_milestoning_file = open(test_milestoning_filename, 'w')
    test_milestoning_file.write(test_milestoning_text)
    test_milestoning_file.close()
    print "test_milestoning_filename:", test_milestoning_filename
    test_model = parse_milestoning_file(test_milestoning_filename)
    self.assertEqual(test_model.num_milestones, 2)
    self.assertEqual(test_model.num_sites, 1)
    self.assertEqual(test_model.sites[0].milestones[0].anchor, '(0.0,0.0,-40.68)')
    self.assertEqual(test_model.sites[0].milestones[1].index, '1')
    self.assertEqual(test_model.temperature, '300')

if __name__ == "__main__": main()
