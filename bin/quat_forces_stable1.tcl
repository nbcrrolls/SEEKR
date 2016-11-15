# quat_forces.tcl

# NAMD tcl script designed to maintain ligand orientation along a particular rotational milestone
# two forces will be imposed maintain ligand rotation along the milestone
#  1. always maintaining a 0.0 dot product from CROSS_QUAT
#  2. if the rotation falls closer to an adjacent milestone's anchor than the starting ANCHOR_QUAT; impose a force opposite this


set outputprefix "SEEKR: "
#set ANCHOR_QUAT "1 0 0 0" ;# the quaternion that represents the center of this rotational milestone
set old_d_theta 0.0;
# NOTE: These must all be in the same hemisphere as ANCHOR_QUAT: that is, all dots with Anchor_quat must be positive
#set NEIGHBOR_QUATS "{0.5 0.5 0.5 0.5}
#                    {0.5 -0.5 -0.5 -0.5}
#                    {0.5 -0.5 0.5 0.5}
#                    {0.5 0.5 -0.5 0.5}
#                    {0.5 0.5 0.5 -0.5}
#                    {0.5 0.5 -0.5 -0.5}
#                    {0.5 -0.5 -0.5 0.5}
#                    {0.5 -0.5 0.5 -0.5}
#                    {0 1 0 0}
#                    {0 0 1 0}
#                    {0 0 0 1}"
#set ANCHOR_QUAT [lindex $NEIGHBOR_QUATS 0]
#set WHICH_NEIGHBOR 1
#set MAIN_NEIGHBOR [lindex $NEIGHBOR_QUATS $WHICH_NEIGHBOR]
set NULL_QUAT "1 0 0 0"
#set EQ_DOT 0.0
#set EQ_ANGLE [expr "2.0 * acos($EQ_DOT)"] ;# 0.0 for tesseract
set PI 3.1415926
set MAX_TORQUE_MAGNITUDE 5.0
#set FORCE_CONSTANT 30.0 ;# in kcal/mol*rad^2
#source /extra/banzai/lvotapka/projects/seekr/old_tcl_functions.tcl
#set LIG_INDECES {1 to 18}
#set REC_INDECES {1 to 18}
# NOTE: need to set recrange
set LA_src "/home/lvotapka/src/vmdplugins/la1.0/la.tcl" ;#NOTE: this will need to be defined by the user
source $LA_src ;# import the tcl linear algebra package

proc assertEqual { value1 value2 } {
  global passes; global failures
  if { $value1 == $value2 } {
    puts " passed"
    incr passes
  } else {
    puts " FAILED: values not equal"
    puts "  value1: $value1"
    puts "  value2: $value2"
    incr failures
    #error "unit test failed."
  }

}

proc assertClose { value1 value2 {tol 0.00005} } {
  global passes; global failures
  set same "True"
  set n [llength $value1]
  for {set i 0} {$i < $n} {incr i} {
    if {[expr "abs([lindex $value1 $i] - [lindex $value2 $i])"] > $tol} {
      set same "False"
    }
  }
  if { $same == "True" } { ;#[expr "abs($value2 - $value1)"] <= $tol
    puts " passed"
    incr passes
  } else {
    puts " FAILED: values not equal"
    puts "  value1: $value1"
    puts "  value2: $value2"
    #error "unit test failed."
    incr failures
  }
}


proc range { start cutoff finish {step 1} } {

  # If "start" and "finish" aren't integers, do nothing:
  if {[string is integer -strict $start] == 0 || [string is integer -strict $finish] == 0} {
    error "range: Range must contain two integers"
  }

  # "Step" has to be an integer too, and
  # no infinite loops that go nowhere are allowed:
  if {$step == 0 || [string is integer -strict $step] == 0} {
    error "range: Step must be an integer other than zero"
  }

  # Does the range include the last number?
  switch $cutoff {
    "to" {set inclu 1}
    "no" {set inclu 0}
    default {
      error "range: Use \"to\" for an inclusive range, or \"no\" for a noninclusive range"
    }
  }

  # Is the range ascending or descending (or neither)?
  set up [expr {$start <= $finish}]

  # If range is descending and step is positive but
  # doesn't have a "+" sign, change step to negative:
  if {$up == 0 && $step > 0 && [string first "+" $start] != 0} {
    set step [expr "$step * -1"]
  }

  # Initialize list variable and identify
  # class of integer range:
  set ranger [list]
  switch "$up $inclu" {
    "1 1" {set op "<=" } ;# Ascending, inclusive range
    "1 0" {set op "<" } ;# Ascending, noninclusive range
    "0 1" {set op ">=" } ;# Descending, inclusive range
    "0 0" {set op ">" } ;# Descending, noninclusive range
  }


  # Generate a list containing the
  # specified range of integers:
  for {set i $start} "\$i $op $finish" {incr i $step} {
    lappend ranger $i
  }
  return $ranger
}

proc parse_selection {selstring} {
  # given a VMD-like selection string, will parse the string into an explicit list of indeces
  if { [catch { ;# attempt to execute code
    set index_list {}
    for {set i 0} {$i < [llength $selstring]} {incr i} {
      if {[lindex $selstring $i] == "to"} {
        foreach item [range [expr "[lindex $selstring $i-1] + 1"] to [expr "[lindex $selstring $i+1]-1"]] {
          lappend index_list $item
        }
        continue
      } else {
        lappend index_list [lindex $selstring $i]
      }
    }
  }]} { ;# then an error occurred
    error "parse_selection in milestoning.tcl: unable to parse selection text: $selstring"
  }
  return $index_list
}

set liglist [parse_selection $LIG_INDECES]
set reclist [parse_selection $REC_INDECES]

set ligatoms {}
foreach atomid $liglist {
  lappend ligatoms $atomid
}
set recatoms {}
foreach atomid $reclist {
  lappend recatoms $atomid
}
#puts "ligatoms: $ligatoms"
foreach atom $ligatoms {
  addatom $atom
}
foreach atom $recatoms {
  addatom $atom
}

proc output { message } { ;# automatically appends caption to stdout from script
  global outputprefix;
  puts "$outputprefix$message"
}

proc cheap_vecdist {vec1 vec2} { ;# distance between the ends of two vectors
  return [veclength [vecsub $vec2 $vec1]]
}

proc cheap_vecnorm {a} { ;# normalize a vector to be of length 1.0
  set a_len [veclength $a]
  if { $a_len == 0.0 } {return "0 0 0"}
  set scalefactor [expr "1/$a_len"]
  return [vecscale $a $scalefactor]
}

proc get_quat_angle { quat1 quat2 } { ;# find the angle between two quaternions ;# old method that used a dot product
  # WARNING: does not use quat_dot, so the quaternions MUST be in the proper hemisphere
  set quatdot [vecdot $quat1 $quat2] ;# get dot product of two quaternions
  return [expr "2.0 * acos($quatdot)"] ;# arccos of the dot product is the angle
}
# NOTE: these two functions give EXACTLY the same result
#proc get_quat_angle { quat1 quat2 } { ;# find the angle between two quaternions
#  # WARNING: does not use quat_dot, so the quaternions MUST be in the proper hemisphere
#  set q_prime [quat_mult $quat2 [quat_inv $quat1]] ;# get dot product of two quaternions
#  set theta [expr "2.0 * acos([lindex $q_prime 0])"] ;# the angle between the two quats
#  return $theta
#}

proc quat_dot {quat1 quat2} {
  # calculates the dot product between two quaternions and makes sure that it returns the proper CLOSEST quat
  set dot1 [vecdot $quat1 $quat2]
  set dot2 [vecdot $quat1 [vecscale -1.0 $quat2]]
  if {$dot1 > $dot2} {
    return $dot1
  } else {
    return $dot2
  }
}

proc quat_conj {quat} {
  # finds the conjugate of a quaternion
  set q0 [lindex $quat 0]
  set q1 [lindex $quat 1]
  set q2 [lindex $quat 2]
  set q3 [lindex $quat 3]
  set result "$q0 -$q1 -$q2 -$q3"
  return $result
}

proc quat_inv {quat} {
  # inverts a quaternion
  set q_sq [vecdot $quat $quat]
  set q_conj [quat_conj $quat]
  set q0 [expr "[lindex $q_conj 0] / $q_sq"]
  set q1 [expr "[lindex $q_conj 1] / $q_sq"]
  set q2 [expr "[lindex $q_conj 2] / $q_sq"]
  set q3 [expr "[lindex $q_conj 3] / $q_sq"]
  return "$q0 $q1 $q2 $q3"
}

proc quat_mult {quat1 quat2} {
  global ANCHOR_QUAT;
  # multiplies two quaternions
  set quat1_0 [lindex $quat1 0]
  set quat2_0 [lindex $quat2 0]
  set quat1_vec [lrange $quat1 1 end]
  set quat2_vec [lrange $quat2 1 end]
  set quatvec_dot [vecdot $quat1_vec $quat2_vec]
  set quatvec_cross [veccross $quat1_vec $quat2_vec]
  set q0 [expr "$quat1_0 * $quat2_0 - $quatvec_dot"]
  set q1 [expr "$quat1_0 * [lindex $quat2_vec 0] + $quat2_0 * [lindex $quat1_vec 0] + [lindex $quatvec_cross 0]"]
  set q2 [expr "$quat1_0 * [lindex $quat2_vec 1] + $quat2_0 * [lindex $quat1_vec 1] + [lindex $quatvec_cross 1]"]
  set q3 [expr "$quat1_0 * [lindex $quat2_vec 2] + $quat2_0 * [lindex $quat1_vec 2] + [lindex $quatvec_cross 2]"]
  set result "$q0 $q1 $q2 $q3"
  return [quat_to_anchor_hemisphere $ANCHOR_QUAT $result]
}

proc quats_to_torque {quat1 quat2 eq_theta force_coeff} {
  global old_d_theta; global MAX_TORQUE_MAGNITUDE
  # given two quaternions, will calculate the torque needed to maintain a constant equilibrium angle for quat1 to reach quat2
  set q_prime [quat_mult [quat_inv $quat2] $quat1] ;# this is the quaternion to make quat1 rotate to quat2
  output "q_prime: $q_prime"
  set axis [lrange $q_prime 1 end] ;# the axis of rotation between the two
  output "axis: $axis"
  if {[veclength $axis] <= 0.00000001} {return {0.0 0.0 0.0}} ;# no torque
  set theta [get_quat_angle $quat1 $quat2] ;# [expr "2.0 * acos([lindex $q_prime 0])"] ;# the angle between the two quats
  output "torque angle: $theta"
  set d_theta [expr "$theta - $eq_theta"] ;# the difference between the quat angle and the equilibrium angle
  output "d_theta: $d_theta"
  if { [expr "abs($d_theta - $old_d_theta)"] > 0.1 } {
    #output "LARGE D_THETA JUMP: [expr {abs($d_theta - $old_d_theta)}]"
    #abort
  }
  set old_d_theta $d_theta
  set torque_magnitude [expr "-2.0 * $force_coeff * $d_theta"] ;# magnitude of the torque
  #if { $torque_magnitude > $MAX_TORQUE_MAGNITUDE } { set torque_magnitude $MAX_TORQUE_MAGNITUDE }
  #if { $torque_magnitude < [expr "-1.0 * $MAX_TORQUE_MAGNITUDE"] } { set torque_magnitude [expr "-1.0 * $MAX_TORQUE_MAGNITUDE"] }
  #output "torque magnitude: $torque_magnitude"
  set torque_unit [cheap_vecnorm $axis] ;# the direction of the torque
  #set torque_unit [cheap_vecnorm [veccross [lrange $quat2 1 end] [lrange $q_prime 1 end] ]]
  #set torque_unit [lrange $quat2 1 end]
  output "torque direction: $torque_unit"
  set torque [ vecscale $torque_unit $torque_magnitude] ;# scale the direction by the magnitude to get the torque vector
  if { [vecdot $torque [vecscale $axis $d_theta"] ] > 0.0 } { ;# then the torque is moving in the same direction as the axis (bad) (NOT NECESSARILY)
    output "TORQUE MOVING IN UNISON TO AXIS"
    output "torque: $torque"
    output "axis: [vecscale $axis $d_theta]"
    #abort
  }
  return $torque
}

proc quat_to_torque {q_prime theta eq_theta force_coeff} {
  # given one quaternion, will calculate the torque needed to maintain a constant equilibrium angle for quat1 to reach quat2
  set axis [lrange $q_prime 1 end] ;# the axis of rotation between the two
  output "axis: $axis"
  if {[veclength $axis] <= 0.00000001} {return {0.0 0.0 0.0}}
  set d_theta [expr "$theta - $eq_theta"] ;# the difference between the quat angle and the equilibrium angle
  set torque_magnitude [expr "-2.0 * $force_coeff * $d_theta"] ;# magnitude of the torque
  set torque_unit [cheap_vecnorm $axis] ;# the direction of the torque
  return [ vecscale $torque_unit $torque_magnitude] ;# scale the direction by the magnitude to get the torque vector
}

# start all the ligand information
proc getCOM { atoms weights } {
  # calculates the center of mass given a list of coordinates (atoms) and weights
  set comsum "0 0 0"
  set totalmass 0
  set n [llength $weights]
  for {set i 0} {$i < $n} {incr i} {
    set mass [lindex $weights $i]; set coord [lindex $atoms $i]
    set tmp [vecscale $mass $coord]
    set comsum [vecadd $comsum $tmp]
    set totalmass [expr "$totalmass + $mass"]
  }
  set com [vecscale [expr "1.0/$totalmass"] $comsum]
  return $com
}

proc getcoords { atoms } {
  loadcoords coords
  set coordlist {}
  foreach atom $atoms {
    lappend coordlist $coords($atom) ;# append the x,y,z coordinates to this
  }
  return $coordlist
}

proc getweights { atoms } {
  loadmasses masses
  set weightlist {}
  foreach atom $atoms {
    lappend weightlist $masses($atom) ;# append the masses to this
  }
  return $weightlist
}

proc confs_to_quat { atoms ref_atoms com refcom weights } {
  # Given a list of atomic coordinates, the reference coordinates, and weights,
  #  this function will return a quaternion that results in the quaternion needed to rotate
  #  the atomic coordinates to the reference coordinates

  # First, find the mean positions of atoms and ref_atoms
  set n [llength $atoms]; # total number of atoms
  if {$n != [llength $ref_atoms]} {error "Error: wrong number of atoms for quaternion comparison" } ;# error check
  # Get relative coordinates by subtracting all positions by their centers of mass, then transforming all atoms to ref_atoms COM
  set relcoord1 {}; set relcoord2 {}
  for {set i 0} {$i < $n} {incr i} {
    lappend relcoord1 [vecsub [lindex $atoms $i] $com]
    lappend relcoord2 [vecsub [lindex $ref_atoms $i] $refcom]
  }
  # Now we need to construct the skew-symmetric matrices
  set M [ make_quat_ls_matrix $relcoord1 $relcoord2 ]
  set q [ trans_principal_eig $M]

  return $q ;# the quaternion to return
}

proc quats_to_anchor_hemisphere { anchor_quat other_quats } {
  # given an anchor quat and a list of other quats, will generate a list where all the other quats
  # have a non-negative dot product with the anchor
  set new_quats {}
  foreach other_quat $other_quats {
    lappend new_quats [quat_to_anchor_hemisphere $anchor_quat $other_quat]
  }
  return $new_quats
}

proc quat_to_anchor_hemisphere { anchor_quat other_quat } {
  # given an anchor quat and a list of other quats, will generate a list where all the other quats
  # have a non-negative dot product with the anchor
  set other_quat_inv [vecscale $other_quat -1.0]
  if {[vecdot $other_quat $anchor_quat] >= [vecdot $other_quat_inv $anchor_quat] } {
    set new_quat $other_quat
  } else {
    set new_quat $other_quat_inv
  }
  return $new_quat
}

output "ANCHOR_QUAT: $ANCHOR_QUAT"
output "MAIN_NEIGHBOR: $MAIN_NEIGHBOR"
#now we need to calculate the CROSS_QUAT based on the projection from the anchor and the main neighbor
set B1 [vecscale [vecdot $ANCHOR_QUAT $MAIN_NEIGHBOR] $ANCHOR_QUAT]
output "B1: $B1"
#set CROSS_QUAT [cheap_vecnorm [vecsub $MAIN_NEIGHBOR $B1]]
set CROSS_QUAT $MAIN_NEIGHBOR
output "CROSS_QUAT: $CROSS_QUAT"
# construct the sampling plane

set $NEIGHBOR_QUATS [quats_to_anchor_hemisphere $ANCHOR_QUAT $NEIGHBOR_QUATS] ;# fix the neighbors so they are all in the same hemisphere as the anchor

set anchor_to_main_neighbor_quat [quat_mult $MAIN_NEIGHBOR [quat_inv $ANCHOR_QUAT]] ;# the normal points out from the plane
output "main_neighbor: $MAIN_NEIGHBOR"
output "anchor_to_main_neighbor_quat: $anchor_to_main_neighbor_quat"
# construct the list that holds all the quaternions that point from the neighbors to the anchor
set neighbor_to_anchor_quats ""
set neighbor_to_anchor_angles ""
set NUM_NEIGHBORS [llength $NEIGHBOR_QUATS]
set max_quat_dot 0.0
foreach neighbor_quat $NEIGHBOR_QUATS {
  lappend neighbor_to_anchor_quats [quat_mult $ANCHOR_QUAT [quat_inv $neighbor_quat]]
  lappend neighbor_to_anchor_angles [get_quat_angle $neighbor_quat $ANCHOR_QUAT ]
  if { [quat_dot $neighbor_quat $ANCHOR_QUAT] == 1.0 } {continue} ;# then this is ourself, so skip
  if { [quat_dot $neighbor_quat $ANCHOR_QUAT] > $max_quat_dot} {
    set max_quat_dot [vecdot $neighbor_quat $ANCHOR_QUAT]
    set closest_quat_angle [expr "acos($max_quat_dot)"]
  }
}
output "neighbor_to_anchor_angles: $neighbor_to_anchor_angles"
output "NUM_NEIGHBORS: $NUM_NEIGHBORS"
output "max_quat_dot: $max_quat_dot"
output "closest_quat_angle: $closest_quat_angle"

# testing vecdot
#output "vecdot {1 0 0 0} {0.5 0.5 0.5 0.5}: [vecdot {1 0 0 0} {0.5 0.5 0.5 0.5}]"
#output "vecdot {0 0.4 0 0} {0.0 2.5 0.0 0.0}: [vecdot {0 0.4 0 0} {0.0 2.5 0.0 0.0}]"
#output "vecdot {3.5 3.5 0.0} {3.5 -3.5 0.0}: [vecdot {3.5 3.5 0.0} {3.5 -3.5 0.0}]"
#abort

# testing veccross
#output "veccross {1 0 0} {0 1 0}: [veccross {1 0 0} {0 1 0}]"
#output "veccross {3 4 5} {-8 7 -6}: [veccross {3 4 5} {-8 7 -6}]"
#output "veccross {0.5 -4.0 -5.3} {-0.002 1.5 2.0}: [veccross {0.5 -4.0 -5.3} {-0.002 1.5 2.0}]"
#abort

#puts "testing function: quats_to_torque"
#assertEqual [quats_to_torque {1 0 0 0} {1 0 0 0} 0.0 100.0] {0.0 0.0 0.0}
#assertEqual [quats_to_torque {1 0 0 0} {0.5 0.5 0.5 -0.5}] { }
#assertClose [quats_to_torque {1 0 0 0} {0.5 0.5 0.5 0.5} [expr "2.0 * $PI / 3.0"] 100.0] {0.0 0.0 0.0}
#assertClose [quats_to_torque {1 0 0 0} {0.5 0.5 -0.5 0.5} 0.0 100.0] "[expr {-400.0 * $PI / (sqrt(3) * 3.0)}] -[expr {-400.0 * $PI / (sqrt(3) * 3.0)}] [expr {-400.0 * $PI / (sqrt(3) * 3.0)}]"

#puts "testing function: quat_to_torque"
#assertEqual [quat_to_torque {1 0 0 0} 0.0 0.0 100.0] {0.0 0.0 0.0}
#assertEqual [quats_to_torque {1 0 0 0} {0.5 0.5 0.5 -0.5}] { }
#assertClose [quat_to_torque {0.5 0.5 0.5 0.5} 4.0 4.0 100.0] {0.0 0.0 0.0}
#assertClose [quat_to_torque {0.5 0.5 -0.5 0.5} [expr "2.0 * $PI / 3.0"] 0.0 100.0] "[expr {-400.0 * $PI / (sqrt(3) * 3.0)}] -[expr {-400.0 * $PI / (sqrt(3) * 3.0)}] [expr {-400.0 * $PI / (sqrt(3) * 3.0)}]"

#puts "testing function: getCOM"
#assertClose [getCOM_fast {3 3 4} {{0.5 0.0 2.5} {3.5 0.0 2.5} {2.0 0.0 0.5}}] {2.0 0.0 1.7}
#assertClose [getCOM_fast {3 3 4} {{0.5 2.5 0.0} {3.5 2.5 0.0} {2.0 0.5 0.0}} ] {2.0 1.7 0.0}

#puts "testing function trans_principal_eig"
#assertClose [trans_principal_eig {{1 3 5 4} {3 5 6 5} {5 6 8 9} {4 5 9 1}}] { -0.3311174067232806  -0.46078481296559143 -0.6687075238512702  -0.48048815453642285}
#assertClose [trans_principal_eig {{1 2 3 4} {2 5 6 7} {3 6 8 9} {4 7 9 10}}] {0.22593827269074562 0.443221861509019 0.5727878807113497  0.6514754961809408}
#assertClose [trans_principal_eig {{10 8 5 2} {8 11 6 4} {5 6 12 9} {2 4 9 13}}] {0.41683205433077336 0.4943681933486974 0.5719781613791227 0.5046702990992139}

#puts "testing function confs_to_quat"
#set atoms    {{1 1 -1} {2 -3 4} {5 -6 7} {-5 7 0} {-5 6 5}}
#set refatoms {{1 1 -1} {2 -3 4} {5 -6 7} {-5 7 0} {-5 6 5}}
#set weights {2.0 2.7 1.0 1.2 0.5}
#set com [getCOM $atoms $weights]
#set refcom [getCOM $refatoms $weights]
#assertClose [confs_to_quat $atoms $refatoms $com $refcom $weights ] {1.0 0.0 0.0 0.0} ;# no rotation
#set atoms {{-2.013885498046875 -2.112818956375122 0.7542771100997925} {4.031986236572266 -0.8651127815246582 2.7267541885375977} {9.1787748336792 -0.6724526882171631 2.0386807918548584} {-8.197371482849121 2.896355390548706 3.8643691539764404} {-4.601040840148926 6.116539001464844 5.506570816040039}}
#set com [getCOM $atoms $weights]
#assertClose [confs_to_quat $atoms $refatoms $com $refcom $weights ] {0.7662610281769211 0.38313051408846055 -0.47891314261057566 -0.19156525704423027} ;# random rotation with same COM
#set atoms {{-6.513885498046875 -8.812818956375121 3.0542771100997923} {-0.4680137634277344 -7.565112781524658 5.0267541885375975} {4.678774833679199 -7.372452688217163 4.338680791854858} {-12.697371482849121 -3.803644609451294 6.16436915397644} {-9.101040840148926 -0.5834609985351564 7.806570816040039}}
#set com [getCOM $atoms $weights]
#assertClose [confs_to_quat $atoms $refatoms $com $refcom $weights ] {0.7662610281769211 0.38313051408846055 -0.47891314261057566 -0.19156525704423027} ;# random rotation with different COM
# maybe later: rotation to a different internal conformation
#abort


proc calcforces { } {
  global ligatoms; global recatoms; global ANCHOR_QUAT; global CROSS_QUAT; global FORCE_CONSTANT; global NULL_QUAT;
  global EQ_ANGLE; global ligrefcoords; global ligweights; global recrefcoords; global recweights; global closest_quat_angle;
  global NEIGHBOR_QUATS; global RECROT; global ligrefquat; global ligquat_addition; global ligrefCOM; global recrefCOM; global recrefquat
  global old_ligquat; global neighbor_to_anchor_quats; global old_d_theta;
  # Initialize step numbers, atom selections, etc.

  #global loadcoords_time; global getcoords_time; global getCOM_time; global confs_to_quat_time; global postprocessing_time; global torque_time

  set stepnum [getstep] ;# step number
  loadcoords coords
  loadmasses masses
  set ligcoords [getcoords $ligatoms]
  set reccoords [getcoords $recatoms]
  if { $stepnum == 0 } { ; # get quaternion of starting rotation
    set ligrefcoords $ligcoords
    set recrefcoords $reccoords
    set ligweights [getweights $ligatoms]
    set recweights [getweights $recatoms]
    set ligrefCOM [getCOM $ligcoords $ligweights]
    set recrefCOM [getCOM $reccoords $recweights]
    set ligrefquat [confs_to_quat $ligrefcoords $ligcoords $ligrefCOM $ligrefCOM $ligweights] ;# should be NULL_QUAT
    set recrefquat [confs_to_quat $recrefcoords $reccoords $recrefCOM $recrefCOM $recweights] ;# should be NULL_QUAT
    set ligquat_addition [quat_mult $ANCHOR_QUAT [quat_inv $ligrefquat]] ;# the quaternion that points from the reference ligand quaternion to ANCHOR_QUAT
    set EQ_ANGLE [get_quat_angle $ligrefquat $CROSS_QUAT]
    set old_ligquat $ligrefquat
    set old_d_theta 0.0
  }

  set ligCOM [getCOM_fast $ligweights $ligcoords ]
  set recCOM [getCOM_fast $recweights $reccoords ]
  output "ligCOM: $ligCOM"
  output "ligrefCOM: $ligrefCOM"
  #  directly find the quaternion that can interpolate the two configurations of the ligand

  # use Quaternion Least Squares to obtain the rotation quaternion
  set ligquat_raw [confs_to_quat $ligcoords $ligrefcoords $ligCOM $ligrefCOM $ligweights]

  # if there is a receptor specified and its rotations are turned 'on'
  if {$RECROT == "True"} {
    set recquat [quat_to_anchor_hemisphere $ANCHOR_QUAT [confs_to_quat $reccoords $recrefcoords $recCOM $recrefCOM $recweights]] ;# get the rotation that recquat has made
  } else { ;# then we are dealing with absolute ligand rotation: it is not rotating with the receptor
    set recquat $NULL_QUAT
  }

  # NOTE: I need to very carefully make sure that this line is correct...
  set ligquat [quat_to_anchor_hemisphere $ANCHOR_QUAT [quat_mult $ligquat_addition [quat_mult $ligquat_raw [quat_inv $recquat]]]] ;# factor out the effect that the receptor has rotated as well as ANCHOR_QUAT
  if { [vecdot $ligquat $old_ligquat] < 0.9} {
    output "LARGE QUATERNION JUMP! [vecdot $ligquat $old_ligquat]"
    output "ligquat: $ligquat"
    output "old_ligquat: $old_ligquat"
    #abort
  }

  # measure distance between current quaternion & anchor quaternion
  set anchor_quat_angle [get_quat_angle $ligquat $ANCHOR_QUAT]
  # impose 1st force: maintaining a desired angle from cross_quat
  set cross_quat_angle [get_quat_angle $ligquat $CROSS_QUAT]
  output "cross quat angle: $cross_quat_angle"
  output "ligquat: $ligquat"
  set U [expr "$FORCE_CONSTANT * ($cross_quat_angle - $EQ_ANGLE) * ($cross_quat_angle - $EQ_ANGLE)"]
  output "EQ_ANGLE: $EQ_ANGLE"
  set lig_torque1 [quats_to_torque $ligquat $CROSS_QUAT $EQ_ANGLE $FORCE_CONSTANT]
  output "lig_torque1: $lig_torque1"
  #if { ( $old_torque1 == 0 ) } { set $old_torque1 $lig_torque1 }
  # impose 2nd force: maintain a minimum quat distance from ANCHOR_QUAT
  set lig_torque2 {0 0 0}
  output "anchor_quat_angle $anchor_quat_angle"
  if { [expr "abs($anchor_quat_angle)"] > $closest_quat_angle } {

    #for {set i 0} {$i < [llength $NEIGHBOR_QUATS]} {incr i} {
    #  set adj_neighbor_quat [lindex $NEIGHBOR_QUATS $i]
    #  output "adj_neighbor_quat: $adj_neighbor_quat"
    #  output "quat_dot ligquat adj_neighbor_quat ${i}: [quat_dot $ligquat $adj_neighbor_quat] "
      #if { [quat_dot $ligquat $ANCHOR_QUAT] < [quat_dot $ligquat $adj_neighbor_quat]} { ;# if ligand rotation is going out of bounds then we need to force the ligand back towards the anchor
    #    set n_to_a [lindex $neighbor_to_anchor_quats $i]
        #set torque_quat [quat_mult $ANCHOR_QUAT [quat_inv $ligquat]] ;# torque toward the anchor
        #set torque_quat $n_to_a ;# torque from the neighbor to the anchor
        #set torque_angle [expr "2.0 * asin( [quat_dot $ligquat $adj_neighbor_quat] - [quat_dot $ligquat $ANCHOR_QUAT])"]
    #    set this_torque2 [quat_to_torque $torque_quat $torque_angle 0.0 $FORCE_CONSTANT]

         set this_torque2 [quats_to_torque $ligquat $ANCHOR_QUAT $closest_quat_angle $FORCE_CONSTANT] ;#$closest_quat_angle $FORCE_CONSTANT]
         output "this_torque2: $this_torque2"
         output "quat_to_anchor: [quat_mult $ANCHOR_QUAT [quat_inv $ligquat]]" ;# torque toward the anchor
         set lig_torque2 [vecadd $lig_torque2 $this_torque2]

    #    set U [expr "$U + ($FORCE_CONSTANT * $torque_angle * $torque_angle * 0.25) "]

        #set sintorque [expr "sin($torque_angle / 2.0)"]
        #set costorque [expr "cos($torque_angle / 2.0)"]
        #set torque_quatprime [cheap_vecnorm [concat $costorque [expr "$sintorque * [lindex $this_torque 0]"] [expr "$sintorque * [lindex $this_torque 1]"] [expr "$sintorque * [lindex $this_torque 2]"]]]


    #    set n_angle [expr "2.0 * acos([lindex $n_to_a 0])"]
    #    set n_axis [lrange $n_to_a 1 end]
    #    set neighbor_to_anchor_axis_angle [vecscale [cheap_vecnorm $n_axis] $n_angle]
    #    set testdot [vecdot $this_torque2 $neighbor_to_anchor_axis_angle]

    #    if { $testdot < 0.0 } {
    #      output "TORQUE OPPOSITE TO ANCHOR: $testdot"
    #      output "neighbor_to_anchor_quat: $n_to_a"
    #      output "torque_quat: $torque_quat"
    #      output "this_torque: $this_torque2"
    #      output "n_angle: $n_angle"
    #      output "n_axis: $n_axis"
    #      output "neighbor_to_anchor_axis_angle: $neighbor_to_anchor_axis_angle"
    #      abort
    #    }
        #output [lindex $neighbor_to_anchor_quats $i]

        #output "U: $U"
        #output "ligquat: $ligquat"
    #    output "torque2_angle: $torque_angle"
        #output "torque2_quat: $torque_quat"
        #output "cross_quat: $CROSS_QUAT"
        #output "eq_angle: $EQ_ANGLE"
        #output "lig_torque1: $lig_torque1"
        #output "lig_torque2: $lig_torque2"
    #  }
    #}

  }

  if { ( [expr "abs([lindex $ligquat 0])"] < 0.01 ) } {
      output "LIGQUAT COLLINEAR! $ligquat"
      abort
    }

  set lig_torque [vecadd $lig_torque1 $lig_torque2]

  # convert torque to a linear force; impose that linear force on all atoms in the ligand
  #output "stepnum: $stepnum  U: $U  lig_torque: $lig_torque"
  if {[veclength $lig_torque] > 0.000001 } { ;# if the torque vector is not tiny
    output "lig_torque: $lig_torque"
    foreach atom $ligatoms {

      set r [vecsub $coords($atom) $ligCOM]
      set unit_torque [cheap_vecnorm $lig_torque]
      #output "unit_torque: $unit_torque"
      #output "r: $r"
      set rho [vecsub $r [vecscale [vecdot $r $unit_torque] $unit_torque]] ;# find the vector rejection in order to get the "arm" length
      set force_magnitude [ expr "[veclength $lig_torque] * [veclength $rho] * $masses($atom)"]
      set force_direction [ cheap_vecnorm [veccross $rho $lig_torque] ]
      set force [vecscale $force_magnitude $force_direction]
      #if {$atom == 0} { output "rho: $rho  force_mag: $force_magnitude  force_direction: $force_direction"}
      addforce $atom $force
    }
  } else {
    set force {0 0 0}
  }

  set old_ligquat $ligquat
  #output "#############################"
}
