#!/usr/bin/python

# NOTE: This script is deprecated

print "WARNING: namd_template is DEPRECATED. Use namd_inputs.py instead."

from adv_template import *
import datetime

class Namd_template(Adv_template):
  def __init__(self, template_filename, params):
    self.template_string = ''.join(open(template_filename, 'r').readlines()) # first load the template file and join (it will be separated by newlines)
    self.params = params
    self.template_string = self.fix_vars()
    rawoutput = self.parse_commands()
    template_string = string.Template(rawoutput) # create the template
    self.output = template_string.safe_substitute(params) # will substitute every indicated $variable in params. Safely in case template file contains extraneous ($$$) dollar signs
    return
    
  def input_gen(self, filename):
    '''generates input file called (filename) and fills the parameters from the dictionary params'''
    out_file = open(filename, 'w') # open output namd input file
    out_file.write(self.output)
    out_file.close()
    return

default_namd_input_params = {
'base':'defaults',
'date':str(datetime.date.today()),

'inpdir':'',
'inpfilename':'dummy',
'outdir':'',
'inpfilename':'dummy',
'firsttimestep':'0',
'timestep':'2.0',

'amber':'yes',
'parmfile':'dummy.prmtop',
'ambercoor':'dummy.inpcrd',
'readexclusions':'yes',
'scnb':'',

'gromacs':'off',
'grotopfile':'',
'grocoorfile':'',


'coordinates':'dummy.pdb',
'structure':'dummy.psf',
'parameters':'CHARMM_parameter_file',
'paratypexplor':'on',
'paratypecharmm':'off',
'velocities':'',
'binvelocities':'$inpname.restart.vel',
'bincoordinates':'$inpname.restart.coor',
'cwd':'',

'outputname':'$outname',
'binaryoutput':'yes',
'restartname':'$outname.restart',
'restartfreq':'500',
'restartsave':'no',
'binaryrestart':'yes',
'dcdfile':'$outname.dcd',
'dcdfreq':'500',
'dcdunitcell':'',
'veldcdfile':'',
'veldcdfreq':'',
'forcedcdfile':'',
'forcedcdfreq':'',

'outputenergies':'500',
'mergecrossterms':'',
'outputmomenta':'',
'outputpressure':'',
'outputtiming':'',
'usegrouppressure':'no',

'extendedsystem':'$inpname.restart.xsc',

'exclude':'scaled1-4',
'_1_4scaling':'0.833333', # NOTE: the '-' has been changed to a '_'. Also, variable may not begin with numerical char
'cutoff':'14',
'switching':'on',
'switchdist':'12',

'pme':'yes',
'pmegridsizex':'',
'pmegridsizey':'',
'pmegridsizez':'',
'pmegridspacing':'',

'rigidbonds':'all',
'rigidtolerance':'0.0005',

'constraints':'',
'consref':'restrain_backbone_ref.pdb',
'conskfile':'restrain_backbone_ref.pdb',
'conskcol':'O',
'constraintscaling':'1.0',

'cellbasisvector1':'',
'cellbasisvector2':'',
'cellbasisvector3':'',
'cellorigin':'',
'xstfile':'$outname.xst'  ,
'xstfreq':'10000',
'wrapwater':'on',
'wrapall':'off',
'wrapnearest':'on',

'minimization':'off',
'minimize':'5000',
'fixedatoms': 'off',
'fixedatomsfile': '',
'fixedatomscol': '',

'nonbondedfreq':'2',

'langevin':'on',
'langevinfile':'',
'langevincol':'O',
'langevintemp':'310',
'langevindamping':'5',
'langevinhydrogen':'no',

'useflexiblecell':'no',
'langevinpiston':'on',
'langevinpistontarget':'1.01325',
'langevinpistonperiod':'100',
'langevinpistondecay':'50',
'langevinpistontemp':'310',

'tclforces':'off',
'tclforcesscript':'',

'colvars':'off',
'colvarsscript':'',

'pairlistdist':'16',
'stepspercycle':'20',
'fullelectfrequency':'20',
'margin':'0.0',
'usesettle':'off'
}
