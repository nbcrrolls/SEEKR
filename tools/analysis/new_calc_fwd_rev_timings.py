#!/usr/bin/python

'''
new_calc_fwd_rev_timings.py
Benjamin Jagger
Amaro Lab 2016

Calculate an estimate of the length of fwd_rev simulation (in ns)
'''
import string, re, glob
import MDAnalysis as mda

OUTPUT_GLOB='*.dcd'
struct_filename= '../../building/holo.parm7'

def extract_timings(output_filenames):
  #timing_lines=[]
  #timings=[]
  #times=[]
  total_time= 0.0
  total_steps=0
  time= 0.0
  for filename in output_filenames:
    #print "Extracting data from", filename
    u= mda.Universe(struct_filename, filename,topology_format='PRMTOP')
    steps= len(u.trajectory)
   # print filename, ' steps', steps




#for line in open(filename,'r'):
#      if re.match(TIMINGS_COMPILE,line):
#        line_split=line.split()
#        time = re.findall('\d.\d+',line_split[3])[0]
#        time=float(time.strip())
#        steps=re.findall('\d+',line_split[1])[0]
#        steps=int(steps.strip())
    #print "steps=" , steps
#        total_time= total_time+time
    total_steps= total_steps+steps
  print "total steps" , total_steps
#  print "total CPU time=" , (total_time* 0.000277778), "hours"
  print "total simulation time (2fs ts)=" , (total_steps*2*1000*0.000001), "nanosecnds"     
 


  return






#for line in open(filename,'r'):
#      if re.match(TIMINGS_COMPILE,line):
#        line_split=line.split()
#        time = re.findall('\d.\d+',line_split[7])
#        times.extend(time)
#    timing_data= [float(x.strip()) for x in times]
#  return timing_data

def main():
  output_filenames= glob.glob(OUTPUT_GLOB)
  timing_data= extract_timings(output_filenames)

if __name__ == "__main__": main()

