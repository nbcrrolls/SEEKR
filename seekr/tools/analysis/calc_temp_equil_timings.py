#!/usr/bin/python

import string, re, glob

TIMINGS_COMPILE= re.compile('TIMING:')
OUTPUT_GLOB='*.out'
def extract_timings(output_filenames):
  #timing_lines=[]
  #timings=[]
  #times=[]
  total_time= 0.0
  for filename in output_filenames:
    print "Extracting data from", filename
    for line in open(filename,'r'):
      if re.match(TIMINGS_COMPILE,line):
        line_split=line.split()
        time = re.findall('\d.\d+',line_split[3])[0]
        time=float(time.strip())
    print "steps=" , steps
    total_time= total_time+time
    total_steps= total_steps+steps
  print "total steps" , total_steps
  print "total CPU time=" , (total_time* 0.000277778), "hours"
  print "total simulation time (2fs ts)=" , (total_steps*2*0.000001), "nanosecnds"
  return

def main():
  output_filenames= glob.glob(OUTPUT_GLOB)
  timing_data= extract_timings(output_filenames)

if __name__ == "__main__": main()

