#!/usr/bin/python

'''
Part of the SEEKR kinetic rate estimation package

contains all the messy information about NAMD min,equil,ensemble, and production simulation parameters

'''
import datetime, os
from adv_template import *
from math import sqrt
import pdb2 as pdb
from pprint import pprint

self_path = os.path.dirname(os.path.realpath(__file__)) # get the path to this script

ambersim_input_template_location = os.path.join(self_path, 'amber_input.template')
colvar_input_template_location = os.path.join(self_path, 'colvar.template')
colvar_script_name = 'colvars.script'

class Ambersim_template(Adv_template):
  def __init__(self, template_filename, params):
    self.template_string = ''.join(open(template_filename, 'r').readlines()) # first load the template file and make it a string separated by newlines
    self.params = params
    print 'template name', template_filename 
    pprint(self.params)
    self.template_string = self.fix_vars()
    rawoutput = self.parse_commands()
    template_string = string.Template(rawoutput) # create the template
    #pprint(template_string)
    self.output = template_string.safe_substitute(params) # will substitute every indicated $variable in params. Safely in case template file contains extraneous ($$$) dollar signs

  def input_gen(self, filename):
    '''generates input file called (filename) and fills the parameters from the dictionary params'''
    out_file = open(filename, 'w') # open output namd input file
    out_file.write(self.output)
    out_file.close()
    return


# Parameters from each dictionary updated in the final NAMD parameters in the order defined below...
default_ambersim_input_params = {
'caption':'default/',
'date':str(datetime.date.today()),

#'inpdir':'',
#'inpfilename':'dummy',
#'outdir':'',
#'outputname':'dummy',
#'firsttimestep':'0',
'dt':'0.002',
'nstlim':'1',
'irest':'0',

#'amber':'yes',
#'parmfile':'dummy.prmtop',
#'ambercoor':'dummy.inpcrd',
#'readexclusions':'yes',
#'scnb':'2.0',
#'exclude':'scaled1-4',
#'_1_4scaling':'0.833333', # NOTE: the '-' has been changed to a '_'. Also, variable may not begin with numerical char
#'watermodel':'',

#'gromacs':'off',
#'grotopfile':'',
#'grocoorfile':'',


#'coordinates':'',
#'structure':'',
#'parameters':'CHARMM_parameter_file',
#'paratypexplor':'on',
#'paratypecharmm':'off',
#'velocities':'',
#'binvelocities':'$inpname.restart.vel',
#'bincoordinates':'$inpname.restart.coor',
#'cwd':'',
'tempi':'298.0',
'temp0':'298.0',
#'watermodel':'',

#'outfilename':'$outname',
#'binaryoutput':'yes',
#'restartname':'$outname.restart',
#'restartfreq':'1000',
#'restartsave':'no',
#'binaryrestart':'yes',
#'dcdfile':'$outname.dcd',
'ntwx':'1000',
#'dcdunitcell':'',
#'veldcdfile':'',
#'veldcdfreq':'',
#'forcedcdfile':'',
#'forcedcdfreq':'',

#'outputenergies':'1000',
#'mergecrossterms':'',
#'outputmomenta':'',
#'outputpressure':'',
#'outputtiming':'',
#'usegrouppressure':'no',

#'extendedsystem':'$inpname.restart.xsc',


'cut':'8',
#'switching':'off',
#'switchdist':'',
#'zeromomentum': 'on',
#'ljcorrection': '',

'ntb':'1',
'ntp':'0',
#'pme':'yes',
#'pmegridsizex':'',
#'pmegridsizey':'',
#'pmegridsizez':'',
#'pmegridspacing':'',

'ntc':'2',
'ntf':'2',
#'rigiditerations':'100',
'tol':'1.0e-8',
'jfastw':'0',
#'usesettle':'on',

#'constraints':'',
#'consref':'restrain_backbone_ref.pdb',
#'conskfile':'restrain_backbone_ref.pdb',
#'conskcol':'O',
#'constraintscaling':'1.0',

#'cellbasisvector1':'',
#'cellbasisvector2':'',
#'cellbasisvector3':'',
#'cellorigin':'',
#'xstfile':'$outname.xst'  ,
#'xstfreq':'10000',
#'wrapwater':'on',
#'wrapall':'off',
#'wrapnearest':'on',
'iwrap':'0',

'imin':'0',
'ntmin':'1',
'maxcyc':'5000',
'ncyc':'1000',
'nmropt': '0',
'ntr':'0',
#'fixedatomsfile': '',
#'fixedatomscol': '',

#'nonbondedfreq':'1',
#'fullelectfrequency':'1',


'ntt':'3',
'gamma_ln':'1.0',
#'langevin':'on',
#'langevinfile':'',
#'langevincol':'O',
#'langevintemp':'300',
#'langevindamping':'5',
#'langevinhydrogen':'no',

#'useflexiblecell':'no',
#'langevinpiston':'on',
#'langevinpistontarget':'1.01325',
#'langevinpistonperiod':'100',
#'langevinpistondecay':'50',
#'langevinpistontemp':'300',
#'useconstantarea':'no',

#'tclforces':'off',
#'tclforcesscript':'',

#'colvars':'off',
#'colvarsscript':'',
#'colvarsconfig':'',


#'pairlistdist':'11',
#'stepspercycle':'20',
#'margin':'0.0',

}

charmm_params = { 
  'amber':'no',
  '_1_4scaling' : '1.0',
  'paraTypeXplor': 'on',
  'cutoff' : '10.0',
  'switchdist' : '10.0',
  'pairlistdist' : '14.0',
}

amber_params = { # based on "Using the Amber force field in NAMD" at http://ambermd.org/namd/namd_amber.html by G. Giambasu and D. Case
  'amber':'yes',
  'readexclusions': 'yes',
  'exclude':'scaled1-4',
  '_1_4scaling' : '0.833333',
  'scnb':'2.0',
  'switching':'off',
  #'switchdist':'12', # if you ever turn on switching, you'll need to define these values
  #'pairlistdist':'11',
  'cutoff':'8',
  'fullelectfrequency' : '1',
  'stepspercycle':'10',
  'ljcorrection':'on',
}

globular_params = {
  'timestep' : '2.0',
  'rigidbonds' : 'all',
  'useflexiblecell':'no',
  'usesettle' : 'on', # since rigidbonds are on, faster than the SHAKE algorithm
  'useconstantarea': 'no',
  'pmegridspacing':'1.0',
}

membrane_params = {
  'timestep' : '2.0',
  'useflexiblecell' : 'yes',
  'usesettle' : 'on',
  'rigidbonds' : 'water',
  'fullelectfrequency' : '1', # frequency (in timesteps) that electrostatics are evaluated
  'useconstantarea': 'no',
  'pmegridspacing':'1.0',
}

min_params = {
  'caption':'SEEKR minimization', # the string at the top of the file
  'imin':'1',
  #'maxcyc':'5000',
  #'pme':'no',
  #'binvelocities':'',
  #'extendedsystem':'',
  #'langevin':'off',
  #'langevinpiston':'off',
  #'timestep':'1.0',
  #'useflexiblecell':'no',
  'iwrap':'0',
  # addition information needed: minimize

}

temperature_equil_params = {
  'caption':'SEEKR temperature equilibration', # the string at the top of the file
  'temp_equil':'True',
  'ntt':3,
  #'tempi':0,
  
  #'pme':'yes',
  #'fullelectfrequency':'1', # frequency (in timesteps) that electrostatics are evaluated
  #'veldcdfile':'$outname.veldcd',
  #'veldcdfreq':'100',
}

ens_equil_params = {
  'caption':'SEEKR ensemble equilibration', # the string at the top of the file
  'irest':'1',
  #'ntp':'0',
  #'barostat':'1',
  #'ntt':'1',
  #'taup':'10.0',
  'ntx':'5',
  #'pme':'yes',
  #'fullelectfrequency':'1', # frequency (in timesteps) that electrostatics are evaluated
  #'veldcdfile':'${outname}vel.dcd',
  #'veldcdfreq':'1000',
}

prod_params = {
  'caption':'SEEKR forward-reverse', # the string at the top of the file
  'temperature':'',
  'pme':'yes',
}

fixed_params = {
  'fixedatoms':'on',
  'fixedatomsfile':'fxd1.pdb',
  'fixedatomscol':'O',
  'zeromomentum':'off', # zeromomentum cannot be on with fixed atoms
}

constraint_params = {
  'constraints' : 'on',
  'consref' : 'constrained1.pdb',
  'conskfile' : 'constrained1.pdb',
  'conskcol' : 'O',
  'zeromomentum':'off', # maybe zeromomentum should be off with harmonic restraints
}

def ensemble_params(ensemble, temp):
  '''populates the ensemble-relevant variables with ensemble information'''
  params = {}
 #if ensemble == 'npt' or ensemble == 'nvt':
  #  params['ntt'] = '1'
    #params['langevintemp'] = temp
    #params['langevindamping'] = '5'
    #params['langevinhydrogen'] = 'no'

  if ensemble == 'npt':
    #params['usegrouppressure'] = 'no'
    params['ntt'] = '3'
    params['ntb']= '2'
    params['ntp']= '1'
    params['gamma_ln']= '1.0'
    params['barostat'] = '1'
    params['pres0'] = '1.01325'
    params['taup'] = '10.0'
    #params['langevinpistonperiod'] = '100'
    #params['langevinpistondecay'] = '50'
    #params['langevinpistontemp'] = temp
  elif ensemble == 'nvt':
    params['ntt'] = '1'
    params['ntb']= '1'
    params['ntp']= '0'
    params['taup']= 10.0
    #params['langevinpiston'] = 'no'
    #params['usegrouppressure'] = 'no'

  elif ensemble == 'nve':
    params['langevin'] = 'off'
    params['langevinpiston'] = 'off'
    params['usegrouppressure'] = 'no'
  else:
    raise Exception, "%s is not a valid ensemble option. Must be 'npt', 'nvt', or 'nve'." % (ensemble,)
  return params



def write_freq_params(freq):
  '''specifies the dcd, xst, ... write frequency'''
  params = {}
  #params['xstfreq'] = freq
  params['ntwx'] = freq
  params['ntwr'] = freq
  params['ntpr'] = freq
  return params

def cell_params(struct, shape="box"):
  '''returns the cellorigin and cellbasisvector parameters for a namd input file'''
  params = {}
  boxdims = pdb.minmax_width(struct) # measures the width of the waterbox
  boxctr =  pdb.center(struct)
  if shape == "box":
    params['cellbasisvector1'] = "%8f 0.0000000 0.0000000" % boxdims[0]
    params['cellbasisvector2'] = "0.0000000 %8f 0.0000000" % boxdims[1] # writes these values to 8 digits
    params['cellbasisvector3'] = "0.0000000 0.0000000 %8f" % boxdims[2]
  elif shape == "oct":
    d = boxdims[0] # assuming that the octahedron is aligned to the x-axis
    params['cellbasisvector1'] = "%8f 0.0000000 0.0000000" % d
    params['cellbasisvector2'] = "%8f %8f 0.0000000" % ((-1.0/3.0)*d, (2.0/3.0)*sqrt(2.0)*d) # writes these values to 8 digits
    params['cellbasisvector3'] = "%8f %8f %8f" % ((-1.0/3.0)*d, (-1.0/3.0)*sqrt(2.0)*d, (-1.0/3.0)*sqrt(6.0)*d)
  else:
    raise Exception, "%s is not a valid ensemble option. Must be 'box' or 'oct'." % (shape,)
   
  params['cellorigin'] = '%8f %8f %8f' % boxctr
  return params
  


def make_input(holo, ff, stage, temperature, write_freq, receptor_type='globular', ensemble='nvt', base='namd input', get_cell=False, fixed=False, constraints=False, settings={}):
  '''creates NAMD input files based on the settings specified:
  arguments:
  ff: can be one of 'amber' or 'charmm'
  stage: can be 'min', 'equil', 'ens', 'prod' # might not actually do this one
  temperatures: a list of temperatures to create simulation files for. Useful for temperature equilibration
  write_freq: the frequency that trajectories, energies, etc. are written
  receptor_type: can be 'membrane' or 'globular'
  ensemble: can be 'nve', 'nvt', 'npt'
  base: a descriptive string added to the input file header
  get_cell: whether the structure is parsed to find waterbox dimensions to populate the cellcenter and cellbasisvector params of the input file
  settings: additional namd parameters that will be substituted into the namd input file
  '''
  params = {}
  params.update(default_ambersim_input_params)
  
  # forcefield
  if ff == 'amber':
    params.update(amber_params)
  elif ff == 'charmm':
    params.update(charmm_params)
  else:
    raise Exception, "%s is not a valid ff option. Must be 'amber' or 'charmm'."


  # receptor_type
  if receptor_type == "globular":
    params.update(globular_params)
  elif receptor_type == "membrane":
    params.update(membrane_params)
  else:
    raise Exception, "%s is not a valid receptor_type option. Must be 'globular', or 'membrane'."
  
  # SEEKR stage
  if stage == 'min':
    params.update(min_params)
  if stage in ['temp_equil', 'equil']:
    params.update(temperature_equil_params)
  if stage == 'ens_equil':
    params.update(ens_equil_params)
  if stage in ['prod','reverse','forward','fwd_rev']:
    params.update(prod_params)

  params.update(write_freq_params(write_freq))
  if get_cell: params.update(cell_params(holo, get_cell)) # if the cell box dimensions needs to be specified, then call it
  if fixed: params.update(fixed_params) # if the cell box dimensions needs to be specified, then call it
  if constraints: params.update(constraint_params) # if the cell box dimensions needs to be specified, then call it

  params['base'] = base
  
  temperature = str(temperature)
  params['temperature'] = temperature
  # ensemble
  params.update(ensemble_params(ensemble,temperature))
  # user-defined
  params.update(settings) # make sure to update with the user-specified parameters last, so they take precedence
  #print "params after ensemble:", params
  ambersim = Ambersim_template(ambersim_input_template_location, params) # generate the namd input class
  return (ambersim, params) # return the input file
  

def read_input(filename):
  ''' reads NAMD input file and returns a dictionary of parameters'''
  paramdict = {}
  namdfile = open(filename,'r')
  for line in namdfile:
    line = line.strip().split()
    if len(line) != 2:
      continue
    key = line[0].lower()

    if key in default_ambersim_input_params.keys(): # then update the value
      paramdict[key] = line[1]
    else:
      pass

  return paramdict

