#!/usr/bin/python

'''
control.py

Controls the parallel processes that run the SEEKR milestoning

Usage:
python control.py command

Commands:
 command [all | anchor #'s] [stage] str_command
 submit [all | anchor #'s] [stage]
 resubmit [all | anchor #'s] [stage]
 check [all | anchor #'s] [stage]
 cancel [all | anchor #'s] [stage]
 modify [all | anchor #'s] [stage] parameter value
 prep [all | anchor #'s] [stage]



Example:
  python control.py submit all min
  python control.py check 6 temp_equil
  python control.py modify all ens_equil numsteps 20000
  python control.py resubmit all ens_equil

'''

import os, unittest, sys, glob, pickle, re, shutil
#import MDAnalysis as mda

try:
  from subprocess import check_output # executes a shell command and returns the STDOUT
  oldpython = False
except ImportError: # then we cannot use this library of subprocess
  oldpython = True

#from sys import argv
from adv_template import Adv_template
from namd_inputs import Namd_template
import namd_inputs
import re
from pprint import pprint
import argparse

namd_inputs_dir = os.path.dirname(os.path.realpath(namd_inputs.__file__))



#print "namd inputs dir:", namd_inputs_dir

#acct = "TG-CHE060073N"
#acct = "TG-MCB140064" # My personal acct number
acct="NONE"
#queue = "normal"
#namd_abrams = "/home1/01624/lvotapka/NAMD_2.9_Source/Linux-x86_64-g++/namd2"
#namd_abrams = "/opt/namd/bin/namd2"
#namd_abrams = "/home1/01624/lvotapka/NAMD_2.9_Source/Linux-x86_64-g++/namd2"
#namd_special = '/home1/03918/tg832177/NAMD_CVS-2015-11-02_Source/Linux-x86_64-g++/namd2' # a special build of NAMD for this purpose
#charm_special = '/home1/03918/tg832177/NAMD_CVS-2015-11-02_Source/Linux-x86_64-g++/charmrun' # required for namd_special above
#mpiexec = '/home1/03918/tg832177/mpiexec'

#program_paths_file='program_paths.pkl'
#program_paths=pickle.load(open(program_paths_file,'rb'))

#namd_special=program_paths['namd_special']
#charm_special=program_paths['charm_special']
#mpiexec=program_paths['mpiexec']

namd_special=os.environ['NAMD_SEEKR']
charm_special=os.environ['CHARM_SEEKR']
mpiexec=os.environ['MPIEXEC']


ALL_COMMANDS = ['command','submit','resubmit','check','cancel','modify','prep']
POSSIBLE_STAGES=['min','temp_equil','ens_equil','fwd_rev',] # the possible stages to search for
DEFAULT_NUM_PROCS = 256 # the default number of processors to assign to a stage job
DEFAULT_TIME_STR = '48:00:00'
PICKLE_NAME = 'cmd.pickle' # the pickle filename: this stores the system object heirarchy in an easily-openable format
TAIL_SIZE = 120 # the number of lines to obtain when tailing a namd output
problem_keys = {'_1_4scaling':'1-4scaling'}

# IDEA: automatically detect the host, or batch processing program, and assign these scripts based on that...

slurm_submit = 'sbatch'
pbs_submit = 'qsub'
submit_cmd = slurm_submit
stampede_submit_template = '''#!/bin/bash
#SBATCH -J  ${job_name}
#SBATCH -o  ${job_name}.%j.out
#SBATCH -e  ${job_name}.%j.err
#SBATCH -p  ${queue}
#SBATCH -n  ${procs}
#SBATCH -t  ${time_str}
#SBATCH -A  ${acct}

module load namd

CUR_PROC=0
echo "starting ${sys_name} ${stage} simulations"

ibrun namd2 ${namd_script} > ${namd_output}

'''

stampede_submit_replica_template = '''#!/bin/bash
#SBATCH -J  ${job_name}
#SBATCH -o  ${job_name}.%j.out
#SBATCH -e  ${job_name}.%j.err
#SBATCH -p  ${queue}
#SBATCH -n  ${procs}
#SBATCH -t  ${time_str}
#SBATCH -A  ${acct}

module load namd
module load fftw3

CUR_PROC=0
echo "starting ${sys_name} ${stage} simulations"

${charm_special} +p${procs} ++scalable-start ++mpiexec ++remote-shell ${mpiexec} ${namd_special} +replicas ${num_replicas} ${namd_script} +stdout ${namd_output}.%d

'''

stampede_submit_generic_template = '''#!/bin/bash
#SBATCH -J  ${job_name}
#SBATCH -o  ${job_name}.%j.out
#SBATCH -e  ${job_name}.%j.err
#SBATCH -p  ${queue}
#SBATCH -n  ${procs}
#SBATCH -t  ${time_str}
#SBATCH -A  ${acct}

echo "starting: ${program}"

ibrun ${program}


'''

gordon_submit_template = '''#!/bin/bash
#PBS -N  ${job_name}
#PBS -o  ${job_name}_${procs_per_node}.out
#PBS -e  ${job_name}_${procs_per_node}.err
#PBS -q  ${queue}
#PBS -l nodes=${nodes}:ppn=${procs_per_node}:
#PBS -l walltime=${time_str}
#PBS -A  ${acct}REPARGS
#PBS -V

module load namd
module load python
module load scipy

cd ${curdir}

echo "starting ${sys_name} ${stage} simulations"

ibrun /opt/namd/2.10b1/bin/namd2 ${namd_script} > ${namd_output}

'''

gordon_submit_generic_template = '''#!/bin/bash
#PBS -N  ${job_name}
#PBS -o  ${job_name}_${procs_per_node}.out
#PBS -e  ${job_name}_${procs_per_node}.err
#PBS -q  ${queue}
#PBS -l nodes=${nodes}:ppn=${procs_per_node}:
#PBS -l walltime=${time_str}
#PBS -A  ${acct}
#PBS -V

module load namd
module load python
module load scipy
module load amber
echo "starting: ${program}"

cd ${curdir}
ibrun ${program}


'''

gordon_submit_serial_template = '''#!/bin/bash
################################################################################
#  Submit script for SDSC's Python-based job bundler for XSEDE/SDSC Gordon.
#  This script shows how to bundle a massive amount of short-running, serial
#  tasks into a single job submission.
#
#  Glenn K. Lockwood, San Diego Supercomputer Center             November 2013
################################################################################
#PBS -N ${job_name}
#PBS -q ${queue}
#PBS -l nodes=${nodes}:ppn=${procs_per_node}:native
#PBS -l walltime=${time_str}
#PBS -v Catalina_maxhops=None,QOS=0
#PBS -o ${job_name}_${procs_per_node}.out
#PBS -e ${job_name}_${procs_per_node}.err
#PBS -A ${acct}
#PBS -V

TASKS=control_file         # the name of your tasks list

cd $PBS_O_WORKDIR

module load python
module load namd
module load scipy
module load amber

pbsdsh /usr/bin/env BUNDLER_MODE=pbsdsh \\
                    PBS_NP=$PBS_NP \\
                    PATH=$PATH \\
                    PYTHONPATH=$PYTHONPATH \\
                    LD_LIBRARY_PATH=$LD_LIBRARY_PATH \\

mpirun_rsh -export \
    -np $PBS_NP \
    -hostfile $PBS_NODEFILE \
    /home/diag/opt/mpi4py/mvapich2/intel/1.3.1/lib/python/mpi4py/bin/python-mpi \
    /home/diag/opt/sdsc-user/bundler/bundler.py $TASKS
'''
# the parameters that are assigned to each SEEKR phase


# Stampede submission parameters
min_params = {
'procs':DEFAULT_NUM_PROCS,'time_str':'12:00:00', 'template':stampede_submit_template, 'acct': acct, 'queue':'normal',
}

temp_equil_params = {
'procs':DEFAULT_NUM_PROCS,'time_str':'02:00:00', 'template':stampede_submit_template, 'acct': acct, 'queue':'normal',
}

ens_equil_params = {
'procs':DEFAULT_NUM_PROCS,'time_str':DEFAULT_TIME_STR, 'template':stampede_submit_template, 'acct': acct, 'queue':'normal',
}

fwd_rev_params = {
'procs':DEFAULT_NUM_PROCS,'time_str':DEFAULT_TIME_STR, 'charm_special': charm_special, 'namd_special':namd_special, 'mpiexec': mpiexec,  'num_replicas':16, 'template':stampede_submit_replica_template, 'acct': acct", 'queue':'normal',
}

submission_template = stampede_submit_template

'''
# Gordon submission parameters
min_params = {
'nodes':1, 'procs_per_node':16,'time_str':'12:00:00', 'template':gordon_submit_template
}

temp_equil_params = {
'nodes':1,'procs_per_node':16,'time_str':'02:00:00', 'template':gordon_submit_template
}

ens_equil_params = {
'nodes':4,'procs_per_node':16,'time_str':'48:00:00', 'template':gordon_submit_template
}

reverse_params = {
'nodes':4,'procs_per_node':16,'time_str':'48:00:00', 'template':gordon_submit_serial_template
}

forward_params = {
'nodes':2,'procs_per_node':16,'time_str':'48:00:00', 'template':gordon_submit_serial_template
}

reverse_prep_params = {
'nodes':1,'procs_per_node':16,'time_str':'24:00:00', 'template':gordon_submit_serial_template
}

forward_prep_params = {
'nodes':1,'procs_per_node':16,'time_str':'6:00:00', 'template':gordon_submit_serial_template
}

submission_template = gordon_submit_template
submit_serial_template = gordon_submit_serial_template
'''


class System():
  '''contains the information concerning the entire system, and how its running'''

  def __init__(self, name, root):
    self.name = name
    self.root = os.path.abspath(root)
    self.anchors=[]

  def find_anchors(self, prefix='anchor'):
    '''searches the system's root directory for all anchors in order to construct the entire tree'''
    print "Pickle not found. Scanning SEEKR root directory for anchors..."
    os.chdir(self.root)
    for anchor_dir in glob.glob('%s*/' % prefix):
      if os.path.isdir(os.path.join(anchor_dir,'md')):
        #print "Found anchor:", anchor_dir
        newanchor = Anchor(anchor_dir)
        newanchor.find_stages(self)
        self.anchors.append(newanchor)
        #print "now returning to directory", self.root
        os.chdir(self.root)

  def get_anchor_by_number(self, number):
    ''' given a numerical value representing the anchor, will return the first anchor that satisfies'''
    for anchor in self.anchors:
      anchor_number = anchor.name.split('_')[1]
      if anchor_number == number: # then we have a match, return this anchor
        return anchor

  def do(self, task, number, stage, defined_params):
    ''' performs a task inside a specific anchor and stage'''
    print "performing a %s" % task
    if task == 'prep': # then remove the control_file that exists in the root directory
      if os.path.isfile(os.path.join(self.root,'control_file')):
        os.remove(os.path.join(self.root, "control_file"))

    if number == 'all':
      for anchor in self.anchors: # loop thru them all
        anchor.do(task,stage, defined_params)
      if task == 'prep': # then this needs to be handled differently
        # SBATCH from this location
        pass
    else:
      # parse by commas and hyphens to generate number_list
      number_list = []
      term_list = number.split(',')
      for term in term_list:
        termsplit = term.split('-')
        if len(termsplit) > 1:
          termsplit = map(str, range(int(termsplit[0]), int(termsplit[1])+1))
        number_list += termsplit
      # now that number_list has been generated
      for num in number_list:
        anchor = self.get_anchor_by_number(num)
        assert anchor != None, "No anchor found to match criteria: %s" % num
        anchor.do(task,stage, defined_params)

class Anchor():
  '''the structure that contains information concerning the anchor'''

  def __init__(self, name ):
    self.name = name
    self.stages = []

  def find_stages(self, system):
    '''finds all stages within this anchor to construct tree object'''
    #print "ls:", os.listdir('.')
    os.chdir(self.name)
    for stage in POSSIBLE_STAGES:
      if not os.path.isdir(os.path.join('md',stage)):
        print "Alert: stage %s does not exist in anchor %s" % (stage, self.name)
        continue
      #print "Found stage:", stage
      newstage = Stage(stage, os.path.abspath(os.path.join('md',stage)), system, self.name, )
      self.stages.append(newstage)

  def unpickle(self, pickle_filename):
    '''given a filename 'pickle_filename', will unpack and reconstruct self tree object'''
    pass # NOTE: not yet implemented

  def do(self, task, stage_name, defined_params): # will print the status of this anchor's stages
    ''' performs a task inside this anchor'''
    print "  anchor %s:" % self.name
    stage_names = stage_name.split(',')
    for stage in self.stages:
      if stage_name == "all" or stage.name in stage_names:
        stage.do(task, defined_params)


class Stage():
  '''represents the process of simulation that needs to be submitted/monitored/etc.'''

  def __init__(self, name,our_dir,system,anchor):
    self.name = name
    self.status = "unsubmitted" # unsubmitted, running, stopped, error, rejected
    self.dir = our_dir
    self.number=1 # number of times has been submitted
    self.system = system
    self.anchor = anchor
    self.procs = DEFAULT_NUM_PROCS
    self.time_str = DEFAULT_TIME_STR



  def submit(self, my_params, my_template):
    # make the submission files
    params = fill_params(self, my_params) #{ 'job_name' : self.name+str(self.number), 'queue':queue,'procs':self.procs, 'time_str':self.time_str, 'acct':acct, 'sys_name':self.system, 'stage':self.name, 'namd_script': self.name+str(self.number)+'.namd', 'namd_output':self.name+str(self.number)+'.out' }
    script_name = '%s%d_%d.submit' % (self.name, self.procs, self.number)
    submit_filename = os.path.join(self.dir, script_name)
    print "      writing submit script to location:", submit_filename

#BRJ 4/4  _0_ should be changed to rotational milestone number
    if self.name == "ens_equil":
      params['job_name']=self.name+"_0_"+str(self.number)
    else:
      params['job_name']=self.name+str(self.number)
#BRJ4/4
    params['curdir']=os.path.join(os.getcwd(), self.dir)
    submit = Adv_template(my_template, params)
    submit_file = open(submit_filename,'w')
    submit_file.write(submit.get_output())
    submit_file.close()
    # submit the files
    cmd = '%s %s' % (submit_cmd, submit_filename)
    submit_stdout = self.command(cmd)
#    print "      executing command:", cmd
    self.status = "submitted"

  def resubmit(self, my_params, my_template):
    '''reads the previous submission file to get the current one's number etc., then resubmits the job'''
    # first, find all previous submission files in this directory
    os.chdir(self.dir)
    prev_out_filename, firsttimestep = self.get_last_timestep()
    prev_num = int(re.findall(r".+(\d+).out", prev_out_filename)[0])
    next_num = prev_num + 1
    prev_prev_num = prev_num - 1
    print "      firsttimestep:", firsttimestep
    self.number = next_num
    new_namd_filename = '%s%d.namd' % (self.name,self.number)
    
    if self.name == 'fwd_rev':
      # just copy the existing fwd_rev namd input file
      old_namd_filename = '%s%d.namd' % (self.name,prev_num)
      shutil.copyfile(old_namd_filename, new_namd_filename)
    else:
      # modify the new NAMD input file
      namd_params_filename='namd_parameters.pkl'
      namd_params=pickle.load(open(namd_params_filename,'rb'))
      newparams = new_namd_output(namd_params, self.name, next_num, prev_num, firsttimestep)
      #prev_namd_filename = '%s%d.namd' % (self.name,prev_num)
      #newparams = new_namd_output(prev_namd_filename, self.name, next_num, prev_num, firsttimestep)
      namd = Namd_template(os.path.join(namd_inputs_dir,namd_inputs.namd_input_template_location),newparams)
      namd.save(new_namd_filename)
      pickle.dump(newparams, open(namd_params_filename, 'wb'))    
    
     # make the actual submission
    self.number = next_num
    self.submit(my_params, my_template)
      
    return 0

  def prep(self, my_template):
    '''prepares the reversal and forward stages for submission'''
    # need to create control file
    if self.name == "bd": # first BD phase
      # we have to have run b_surface
      cmd = "cd %s; python extract_bd_frames.py %s\n" % (self.dir, traj_dir) # NOTE: traj_dir doesn't yet exist

      self.status = "prepped"

  def modify(self, param, newval):
    ''' runs through .namd files selects most recent and modifies specified param to newval'''
    os.chdir(self.dir)
    param= param.lower()
    #prev_num = int(re.findall(r".+(\d+).out", prev_out_filename)[0])
    # modify the new NAMD input file
    namd_params_filename='namd_parameters.pkl'
    params=pickle.load(open(namd_params_filename,'rb'))

    if param in params:
      if params[param] != "":
        namd_filename, namd_filenumber = self.get_last_input()
        namd_file= '%s.namd' %(namd_filename)
        print "      NAMD input file modified:" , namd_file 
        params[param] = newval
        namd = Namd_template(os.path.join(namd_inputs_dir,namd_inputs.namd_input_template_location),params)
        namd.save(namd_file)
      else:
        print "      Parameter not used in this stage! No modification made" 
    else:
      print "      Invalid parameter! No modification made"
    #pprint(params)
    pickle.dump(params, open(namd_params_filename, 'wb'))
    return 0

  def command(self, cmd, indir=True):
    '''will execute a specific command within the stage's directory.
    Returns the standard output of that command.'''
    curdir = os.getcwd() # get the directory we are in
    if indir: os.chdir(self.dir) # move to the stage's directory
    if not oldpython:
      std_out = check_output(cmd, shell=True) # execute the command
    else:
      std_out = ''
    #print "command executed with return code:", result # alert the user of the result of the command
    os.chdir(curdir) # return to the original directory
    return std_out

  def get_last_timestep(self, file_index=-1):
    '''finds the last timestep executed for this stage'''
    os.chdir(self.dir)
    output_files = glob.glob('%s*.out*' % self.name) # use glob pattern matching to get all the output files

    output_files = [outfile for outfile in output_files if re.match("%s\d+\.out.*" % self.name, outfile)] # further filter the files from the glob output
    # now the files must be sorted to get the latest one
    #print "output files:", output_files
    if not output_files:
      print "      Alert: no namd output files found."
      return '', "0"
    try:
      prev_out_filename = sorted(output_files, key=lambda item: int(re.findall(r".(\d+)\.out.*",item)[0]))[file_index]
    except IndexError:
      print "IndexError namd output files found, but a problem occurred parsing the output stage numbers."
      print "output_files:", output_files
      print "regular expression: %s not matched." % r".(\d+)\.out.*"
      exit()

    #prev_num = int(re.split(r"_|\.",prev_output)[1])
    #next_num = prev_num + 1
    #prev_prev_num = prev_num - 1
    # now we need to read the previous namd output file in order to get the last step number
    #prev_out_filename = self.name+str(prev_num)+'.out'
    if self.name == "fwd_rev": return prev_out_filename, "0"
    
    if not os.path.exists(prev_out_filename):
      self.status
      print "      Alert: namd output file not found:", prev_out_filename
      return '', "0"
    prev_out_file = open(prev_out_filename,'r')
    tail_list = tail(prev_out_file,TAIL_SIZE)[0]
    tail_list.reverse()
    lastline = None
    for line in tail_list:
      lastline = re.findall("WRITING VELOCITIES TO RESTART FILE AT STEP (\d+)", line)
      if lastline:
        break
    if not lastline:
      print "      trying earlier output file in case last output was an error..."
      return self.get_last_timestep(-2)
    return prev_out_filename, lastline[0]

#BRJ 2/4 4/4

  def get_last_input(self, file_index=-1):
    os.chdir(self.dir)
    input_files = glob.glob('%s*.namd' % self.name) # use glob pattern matching to get all the namd input files
    
    if self.name == "ens_equil":     
      input_files = [infile for infile in input_files if re.match("%s_\d+_\d+\.namd" % self.name, infile)] # further filter the files from the glob output
    elif self.name == "fwd_rev":
      input_files = [infile for infile in input_files if re.match("%s\d+\.namd" % self.name, infile)]
    else:
      print "This stage cannot be modified" 
   # now the files must be sorted to get the latest one
    if not input_files:
      print "      Alert: no namd input files found."
      return '', "0"
    try:
      prev_in_filename = sorted(input_files, key=lambda item: int(re.findall(r".(\d+)\.namd",item)[0]))[file_index]
    except IndexError:
      print "IndexError detected which finding previous NAMD input files"
      exit()
    if not os.path.exists(prev_in_filename):
      self.status
      print "      Alert: namd input file not found:", prev_in_filename
      return '', "0"
    prev_in_num = int(re.findall(r".+(\d+).namd", prev_in_filename)[0])
    line= prev_in_filename.split('.')
    prev_in_filename= line[0]
    return (prev_in_filename, prev_in_num) 


#end BRJ 2/4


  def do(self, task, defined_params):
    ''' performs a task inside this stage'''
    print "    %s" % self.name
    my_template = submission_template
    my_params = {}
    if self.name == 'min':
      my_params = min_params
    if self.name == 'temp_equil':
      my_params = temp_equil_params
    if self.name == 'ens_equil':
      my_params = ens_equil_params
    if self.name == 'fwd_rev':
      my_params = fwd_rev_params
      my_template = my_params['template']

    if defined_params['time']:
      my_params['time_str'] = defined_params['time']

    if defined_params['procs']:
      my_params['procs'] = int(defined_params['procs'])
      self.procs = my_params['procs']

    if defined_params['acct']:
      my_params['acct'] = defined_params['acct']

    if defined_params['queue']:
      my_params['queue'] = defined_params['queue']

    if defined_params['num_replicas']:
      my_params['num_replicas'] = defined_params['num_replicas']

    if defined_params['procs_per_node']:
      my_params['procs_per_node'] = defined_params['procs_per_node']

    if task in ["check", "status"]:
      print "      last timestep: %s" % self.get_last_timestep()[1]
      print "      status: %s" % self.status

    elif task == "submit":
      print "      now submitting job in stage: %s" % self.name
      self.submit(my_params, my_template)

    elif task == "resubmit":
      print "      now resubmitting job in stage: %s" % self.name
      #if self.name == 'fwd_rev':
      #  self.submit(my_params, my_template) # only a submit is necessary
      #else:
      self.resubmit(my_params, my_template)

    elif task == "prep": # NOT used in this version of SEEKR
      print "      Alert: unable to prep job for stage: %s" % self.name

    elif task.startswith("command"):
      cmd = ' '.join(task.split()[1:])
      print "      now executing: '%s' in stage: %s" % (cmd, self.name)
      submit_stdout = self.command(cmd)
      print "      result:", submit_stdout

    elif task == "cancel":
      #print "      now cancelling job launched from stage: %s" % self.name
      print "      Alert: 'cancel' command not yet implemented."

    elif task.startswith("modify"):
      param = task.split()[1]
      newval = task.split()[2]
      print "      now modifying parameter: '%s' in stage: %s to value: %s" % (param, self.name, newval) # NOTE: this needs to be fixed
      self.modify(param,newval)

    else:
      print "      Alert: command %s not recognized" % task


def tail(f, n, offset=None):
  """Reads a n lines from f with an offset of offset lines.  The return
  value is a tuple in the form ``(lines, has_more)`` where `has_more` is
  an indicator that is `True` if there are more lines in the file.
  """
  avg_line_length = 74
  to_read = n + (offset or 0)

  while 1:
    try:
      f.seek(-(avg_line_length * to_read), 2)
    except IOError:
      # woops.  apparently file is smaller than what we want
      # to step back, go to the beginning instead
      f.seek(0)
    pos = f.tell()
    lines = f.read().splitlines()
    if len(lines) >= to_read or pos == 0:
      return lines[-to_read:offset and -offset or None], \
           len(lines) > to_read or pos > 0
    avg_line_length *= 1.3

def new_namd_output(srcparams, stage, next_num, prev_num, firsttimestep):
  '''using an existing namd input file, will create a new one with updated characteristics'''
  #srcparams = namd_inputs.read_input(srcfile)
 # srcparams=pickle.load(open(srcfile,'r+b'))
  srcparams['firsttimestep']=firsttimestep
  srcparams['inpdir']=''
  srcparams['inpfilename']='%s%d' % (stage, prev_num)
  srcparams['outdir']=''
  srcparams['outfilename']='%s%d' % (stage, next_num)
  #pickle.dump(namd_params, (open(srcfile, 'wb'))
  #for problem_key in problem_keys:
  #  if problem_key in srcparams.keys():
  #    srcparams[problem_key] = problem_keys[problem_key]
  #print "srcparams:"
  #pprint(srcparams)
  return srcparams


'''def extract_from_ens(directory):
  pos_dcds = glob.glob("../ens_equil/ens_equil*[0-9].dcd") # get all coordinate dcds from the ens_equil folder
  vel_dcds = glob.glob("../ens_equil/ens_equil*[0-9]*vel.dcd") # get all vel dcds from the ens_equil folder
  misnamed_veldcds = glob.glob("../ens_equil/ens_equil*.veldcd") # if they are misnamed with the extension '.veldcd'
  print "misnamed veldcds:", misnamed_veldcds
  for i in misnamed_veldcds: # must rename all .veldcd to vel.dcd
    newname = re.sub('.veldcd','vel.dcd',i)
    print "oldname: ",i," newname: ", newname
    os.rename('../ens_equil/'+i,'../ens_equil/'+newname) # rename the file in the filesystem
    vel_dcds.append(newname) # append the correct name to vel_dcds list

  pos_dcds.sort()
  vel_dcds.sort()
  print "pos_dcds:", pos_dcds
  print "vel_dcds:", vel_dcds

pos_universe = mda.Universe("../building/holo.prmtop", pos_dcds, format='dcd')
vel_universe = mda.Universe("../building/holo.prmtop", vel_dcds, format='dcd')
pos_all = pos_universe.selectAtoms("all") # select all the atoms in position file
vel_all = vel_universe.selectAtoms("all")
pos_traj_len = pos_universe.trajectory.numframes
vel_traj_len = vel_universe.trajectory.numframes
if pos_traj_len != vel_traj_len:
  print "WARNING: extraction script: position trajectory and velocity trajectory of unequal lengths"
  # raise Exception "position trajectory and velocity trajectory of unequal lengths"
frame_list = range(30000,100001,20) # we should have a list of length num_frames now, spread uniformly across the
for frame_index in frame_list:
  print "now writing frame", frame_index
  pos_universe.trajectory[frame_index] # set the frame
  vel_universe.trajectory[frame_index] # set the frame
  pos_all.write(os.path.join("./","pos%i.pdb" % frame_index)) # write the pdb file of positions at this frame
  vel_all.write(os.path.join("./","vel%i.pdb" % frame_index)) # write the pdb file of velocities at this frame
'''

def fill_params(stage, starting_params):
  # I should code this function to read input from the user more effectively, using ArgParse
  d = {}
  d.update(starting_params)
  #d['acct'] = acct
  #d['queue'] = queue
  d['sys_name'] = stage.system.name
  d['stage'] = stage.name
#BRJ 4/4  _0_ should be changed to rotational milestone number
  if stage.name== "ens_equil":
    d['namd_script'] = stage.name+'_0_'+str(stage.number)+'.namd'
  else:
    d['namd_script'] = stage.name+str(stage.number)+'.namd'
  if stage.name== "ens_equil":
    d['namd_output'] = stage.name+'_0_'+str(stage.number)+'.out'
  else:
    d['namd_output'] = stage.name+str(stage.number)+'.out'
#BRJ 4/4
  if 'procs' in starting_params.keys():
    stage.procs = starting_params['procs']
  else:
    stage.procs = starting_params['procs_per_node'] * starting_params['nodes']
  return d

if __name__ == "__main__":
  # first extract the commands

  parser = argparse.ArgumentParser(description="The SEEKR program that is designed to make the numerous and parallelized MD phases to be easier to submit, manage, prepare, and control on a supercomputer or cluster. Examples: 'python control.py submit all min',  'python control.py check 6 temp_equil',  'python control.py modify all ens_equil numsteps 20000', 'python control.py resubmit all ens_equil'")
  # positional arguments
  command_choices = ["submit", "resubmit", "command", "check", "status", "cancel","modify", "prep" ]
  parser.add_argument('command', metavar='COMMAND', type=str, choices=command_choices, help="the command to execute. Options include: %s" % (', '.join(command_choices)))
  parser.add_argument('anchors', metavar='ANCHORS', type=str, help="the anchors to execute the command for. May be the keyword 'all', a single integer, or comma-separated list of integers, or dash to indicate a number range. Examples: 'all', '3', '2,5,7', '1-4', '1,3-6,8'.")
  parser.add_argument('stage', metavar='STAGE', type=str, help="Which stage of the anchor(s) to apply the command to. At this time, relevant stages are 'ens_equil' and 'fwd_rev'.")
#BRJ 2/1
 
  parser.add_argument('param', metavar='PARAM', type=str, nargs='?', help="the parameter that will be changed in the .namd file")
  parser.add_argument('newval', metavar='NEWVAL', type=str,nargs='?', help="the value that the parameter will be changed to")   

  # optional arguments
  parser.add_argument('-t', '--time', metavar="TIME", dest='time', default='', type=str, help='The time to run the simulations.')
  parser.add_argument('-p', '--procs', metavar="PROCS", dest='procs', default='', type=str, help='The number of processors to request.')
  parser.add_argument('-a', '--acct', metavar="ACCT", dest='acct', default='', type=str, help='The account number to request with.')
  parser.add_argument('-q', '--queue', metavar="QUEUE", dest='queue', default='', type=str, help='The queue to submit to.')
  parser.add_argument('-nr', '--num_replicas', metavar="NUM_REPLICAS", dest='num_replicas', default='', type=str, help='The number of replicas to use for the calculation (if applicable, such as in the fwd_rev phase).')
  parser.add_argument('-ppn', '--procs_per_node', metavar="PROCS_PER_NODE", dest='procs_per_node', default='', type=str, help='If needed by the supercomputer, the number of processors per node to use.')

  sys_root = os.path.abspath('.')
  sys_name = os.path.basename(sys_root)
  system = System(sys_name, sys_root)

  args = parser.parse_args() # parse all the arguments
  args = vars(args)
  command = args['command']
  anchors = args['anchors']
  stage = args['stage']
  param= args['param']
  newval= args['newval']
  timestr = args['time']
  numprocs = args['procs']
  acct = args['acct']
  queue = args['queue']
  num_replicas = args['num_replicas']
  ppn = args['procs_per_node']

  # try to unpack the filetree structure from a pickle file
  task = command
  if command in ["command","modify"]:
    task= command + ' ' + param +' ' + newval

  if os.path.isfile(PICKLE_NAME):
    pickleread = open(PICKLE_NAME,'rb') # open the pickle file
    system = pickle.load(pickleread) # read the pickle
    pickleread.close()

  else: # if the pickle doesn't exist, then create a new one
    system.find_anchors() # since the pickle doesn't exist, then create a new tree
  system.do(task, anchors, stage, args) # now execute whatever command we're given

  # save the pickle
  os.chdir(sys_root)

  '''
  assert len(argv) >= 2, "script requires an argument."
  command = argv[1] # the command
  if command not in ALL_COMMANDS:
    print "ALERT: command unknown: %s. Aborting..." % command
    exit()

  sys_root = os.path.abspath('.')
  sys_name = os.path.basename(sys_root)
  system = System(sys_name, sys_root)

  if command in ["submit", "resubmit", "command", "check", "status", "cancel","modify", "prep" ]:
    assert len(argv) > 3, "command %s requires additional arguments." % command
    anchors = argv[2]
    stage = argv[3]

    task = command
    if command in ["command","modify"]:
      cmd_str = ' '.join(argv[4:]) # concatenate the command to execute
      task = command + ' ' + cmd_str
    #if anchors == "all": # then its being applied to all anchors
      #for anchor in



    # try to unpack the filetree structure from a pickle file
    if os.path.isfile(PICKLE_NAME):
      pickleread = open(PICKLE_NAME,'rb') # open the pickle file
      system = pickle.load(pickleread) # read the pickle
      pickleread.close()

    else: # if the pickle doesn't exist, then create a new one
      system.find_anchors() # since the pickle doesn't exist, then create a new tree
    system.do(task, anchors, stage) # now execute whatever command we're given

    # save the pickle
    os.chdir(sys_root)
    #print "saving the pickle in:", os.path.join(sys_root, PICKLE_NAME)
    #picklewrite = open(PICKLE_NAME,'wb') # open the pickle file for writing
    #pickle.dump(system, picklewrite) # write the pickle
    #picklewrite.close()

  elif command == "help":
    print __doc__ # display help information

  elif command == "repickle":
    system.find_anchors() # recreating the pickle
  '''
