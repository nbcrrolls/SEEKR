# Test script that checks the outputs of all the functions in 'milestoning.tcl' to make sure they are producing the expected results



# a list of dummy variables that need to be defined in order for the functions to be tested
set LIGRANGE {1 to 4 }
set SCRIPT_INTERVAL 1
set LIGROT False
set ABORT_ON_CROSSING False
set PHASE forward
set CARE_ABOUT_SELF False
set MILESTONE_LIST {}
set whoami 0
set GRID_EDGE_RAD 1.0
set RECRANGE {5 to 8}
set RECPA1_LIST {}
set RECPA3_LIST {}
#set FHPD_FILE "../test/fhpd_test.txt"
proc addgroup {arg} {} ;# dummy function

# a list of all the functions to test
set proclist { 
easy_principalaxes 
range 
parse_selection
cheap_vecdist 
cheap_vecnorm 
cheap_vecdot 
mat_mat_mult 
mat_vec_mult 
translate_mat 
cheap_trans 
cheap_vecbackscale 
check_spherical_milestone 
check_planar_milestone
calc_ref_movements
calcforces
}
#add_id_to_FHPD 
#orient

source ../bin/milestoning.tcl

proc Assert_equal { arg1 arg2 {func "not specified"}} {
  if {$arg1 != $arg2} {error "Assert_equal: in function $func, arguments not equal: $arg1, $arg2"} else {return "passed"}
}

#proc test_add_id_to_FHPD {} {
#  global FHPD_FILE
#  file delete $FHPD_FILE
#  add_id_to_FHPD 42
#  add_id_to_FHPD 63
#  add_id_to_FHPD 0
#  set testfile [open $FHPD_FILE r]
#  set lines [read $testfile]
#  set linelist [split $lines]
#  close $testfile
#  Assert_equal $linelist {42 63 0 {}} add_id_to_FHPD
#  return passed
#}

#proc test_orient {} {
#  set vec1 {2 3 4}
#  set vec2 {-6 3 -1}
#  set mat1 [orient {0 0 0} $vec1 $vec2]
#  set finalvec1 [vecnorm [ lrange [mat_vec_mult $mat1 "$vec1 1"] 0 2 ]]
#  set finalvec2 [vecnorm $vec2]
#  if {[vecdot $finalvec1 $finalvec2] < 0.99999} {error "'orient' function not working properly: test failed. vec1: $finalvec1, vec2: $finalvec2"}
#  return passed
#}

proc test_easy_principalaxes {} {
  puts "ALERT: this function has no unit tests."
}

proc test_range {} {
  Assert_equal [range 5 to 9] {5 6 7 8 9}
  return passed
}

proc test_parse_selection {} {
  Assert_equal [parse_selection "3 5 to 9 14 16 to 18"] {3 5 6 7 8 9 14 16 17 18}
}

proc test_cheap_vecdist {} {
  set testvec1 {4 0 0}
  set testvec2 {-8 -5 0}
  Assert_equal [cheap_vecdist $testvec1 $testvec2] 13.0 cheap_vecdist
  return passed
}

proc test_cheap_vecnorm {} {
  puts "ALERT: this function has no unit tests."
}

proc test_cheap_vecdot {} {
  puts "ALERT: this function has no unit tests."
}

proc test_mat_mat_mult {} {
  puts "ALERT: this function has no unit tests."
}

proc test_mat_vec_mult {} {
  puts "ALERT: this function has no unit tests."
}

proc test_translate_mat {} {
  puts "ALERT: this function has no unit tests."
}

proc test_cheap_trans {} {
  puts "ALERT: this function has no unit tests."
}

proc test_cheap_vecbackscale {} {
  puts "ALERT: this function has no unit tests."
}

proc test_check_spherical_milestone {} {
  Assert_equal [check_spherical_milestone {1 2 3} {4 5 6} 10.0] inside
  Assert_equal [check_spherical_milestone {1 2 3} {4 5 -60} 10.0] outside
  Assert_equal [check_spherical_milestone {-1 -2 -3} {-9 -10 -11} 10.0] outside
  return passed
}

proc test_check_planar_milestone {} {
  Assert_equal [check_planar_milestone {3 4 5} {4 5 6} {-1 0 0}] positive
  Assert_equal [check_planar_milestone {0 10 3} {-9 10 3} {1 1 1}] positive
  Assert_equal [check_planar_milestone {1 2 3} {0 0 0} {-1 -1 0}] negative
  return passed
}

proc test_calc_ref_movements {} {
  puts "ALERT: this function has no unit tests."
}

proc test_calcforces {} {
  puts "ALERT: this function has no unit tests."
}

proc main {} {
  global proclist
  foreach func $proclist {
    puts "Now testing function: $func"
    set cmd "test_$func"
    set result [$cmd]
    puts "$result"
  }
}

main
