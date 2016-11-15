# quat_forces.tcl

# NAMD tcl script designed to maintain ligand orientation along a particular rotational milestone
# two forces will be imposed maintain ligand rotation along the milestone
#  1. always maintaining an equal dot product from CROSS_QUAT to ANTICROSS_QUAT
#  2. if the rotation falls closer to an adjacent milestone's anchor than the starting ANCHOR_QUAT; impose a force opposite this


set outputprefix "SEEKR: "
set old_d_theta 0.0
set NULL_QUAT "1 0 0 0"
set PI 3.1415926
set MAX_TORQUE_MAGNITUDE 5.0
#set LA_src "/home/lvotapka/src/vmdplugins/la1.0/la.tcl" ;#NOTE: this will need to be defined by the user
#source $LA_src ;# import the tcl linear algebra package

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
  if {$quatdot > 1.0} {set quatdot 1.0}
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

proc get_S_matrix_transpose { q } {
  # find the transpose of the S-matrix: a matrix that contains quaternion information
  set q0 [lindex $q 0]; set q1 [lindex $q 1]; set q2 [lindex $q 2]; set q3 [lindex $q 3]
  set q0_neg [expr "-1.0 * $q0"]; set q1_neg [expr "-1.0 * $q1"]; set q2_neg [expr "-1.0 * $q2"]; set q3_neg [expr "-1.0 * $q3"];
  set S_t "{ $q0 $q1 $q2 $q3} { $q1_neg $q0 $q3 $q2_neg } { $q2_neg $q3_neg $q0 $q1 } { $q3_neg $q2 $q1_neg $q0 } "
  return $S_t
}

proc get_quat_potential_gradient { q force_coeff {quat_indeces {0 1 2 3}} {coeffs {1.0 1.0 1.0 1.0}} } {
  # get the gradient of the quaternion potential energy function
  # U(q) = force_coeff * (coeff0*q0 + coeff1*q1 + coeff2*q2 + coeff3*q3) ** 2
  set quat_sum 0.0
  for {set i 0} {$i < [llength $quat_indeces]} {incr i} {
    set quat_index [lindex $quat_indeces $i]
    set coeff [lindex $coeffs $i]
    set quat_sum [expr "$quat_sum + ($coeff * [lindex $q $quat_index])"]
    lappend grad_U []
  }
  output "quat_sum: $quat_sum"
  set grad_U {}
  set dU_dqi [expr "2.0 * $force_coeff * ($quat_sum)"]
  for {set i 0} {$i < [llength $quat_indeces]} {incr i} {
    lappend grad_U [expr "[lindex $coeffs $i] * $dU_dqi"]
  }
  #output "grad_U: $grad_U"
  return $grad_U
}

proc get_quat_potential_gradient_angle { q eq_theta force_coeff } {
  # get the gradient of the quaternion potential energy function
  # U(q) = force_coeff * (2.0*acos(q0) - eq_theta) ** 2
  set q0 [lindex $q 0]
  set angle_diff [expr "(2.0*acos($q0) - $eq_theta)"]
  if { [expr "abs($angle_diff)"] < 0.0000001} { return "0.0 0.0 0.0 0.0" }
  set dU_dq0 [expr "4.0 * $force_coeff * $angle_diff / sqrt(1 - $q0 * $q0)"]
  set grad_U "$dU_dq0 0.0 0.0 0.0"
  return $grad_U
}

proc get_quat_potential_gradient_2 { q1 q2 q3 force_coeff } {
  # get the gradient of the quaternion potential energy function
  # U(q1, q2, q3) = force_coeff * (2.0 * acos(dot(q1, q2)) - 2.0 * acos(dot(q1, q3)))**2
  set gradU {}
  set prefix_q1_q2 [expr "8.0 * $force_coeff * (acos([quat_dot $q1 $q2]) - acos([quat_dot $q1 $q3])) / sqrt(1 - [quat_dot $q1 $q2]*[quat_dot $q1 $q2])"]
  set prefix_q1_q3 [expr "8.0 * $force_coeff * (acos([quat_dot $q1 $q2]) - acos([quat_dot $q1 $q3])) / sqrt(1 - [quat_dot $q1 $q3]*[quat_dot $q1 $q3])"]
  for { set i 0 } {$i < 4} {incr i} {
    set q2i [lindex $q2 $i]; set q3i [lindex $q3 $i]
    lappend gradU [expr "$prefix_q1_q2 * $q2i - $prefix_q1_q3 * $q3i"]

  }
  return $gradU
}

proc get_torque_from_quat_potential {q_prime force_coeff {quat_indeces {0 1 2 3} } { coeffs {1 1 1 1} } } {
  # given a quaternion, will calculate the torque imposed by a harmonic potential on a sum of the quaternion multiplied by coeffs
  # torque = -0.5 * S^T(q_prime) * grad(potential(q_prime)) + torque_internal
  set S_t [get_S_matrix_transpose $q_prime]
  set quat_grad [get_quat_potential_gradient $q_prime $force_coeff $quat_indeces $coeffs]
  set t4 [vecscale [vectrans $S_t $quat_grad] 0.5 ] ;# the torque 4-vector
  set torque [lrange $t4 1 end] ;# get the torque 3-vec from the torque 4-vec
  return $torque
}

proc get_torque_from_quat_potential_angle {q_prime eq_theta force_coeff } {
  # given a quaternion, will calculate the torque imposed by a harmonic potential imposed on the angle
  # torque = -0.5 * S^T(q_prime) * grad(potential(q_prime)) + torque_internal
  set S_t [get_S_matrix_transpose $q_prime]
  set quat_grad [get_quat_potential_gradient_angle $q_prime $eq_theta $force_coeff]
  set t4 [vecscale [vectrans $S_t $quat_grad] -0.5 ] ;# the torque 4-vector
  set torque [lrange $t4 1 end] ;# get the torque 3-vec from the torque 4-vec
  return $torque
}

proc get_torque_from_quat_potential_2 {q1 q2 q3 force_coeff} {
  # given a quaternion, will calculate the torque imposed by a harmonic potential imposed on the angle
  # torque = -0.5 * S^T(q_prime) * grad(potential(q_prime)) + torque_internal
  set S_t [get_S_matrix_transpose $q1]
  set quat_grad [get_quat_potential_gradient_2 $q1 $q2 $q3 $force_coeff] ;# gradient of potential energy function
  set t4 [vecscale [vectrans $S_t $quat_grad] -0.5 ] ;# the torque 4-vector
  set torque [lrange $t4 1 end] ;# get the torque 3-vec from the torque 4-vec
  return $torque
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

proc get_moment_of_inertia { weights atoms pivot_vec COM } {
  # calculates the moment of inertia for a molecule
  set n [llength weights]
  set moment 0.0
  for {set i 0} {$i < $n} {incr i} {
    set mass [lindex $weights $i]
    set atom [lindex $atoms $i]
    set r [vecsub $atom $COM]
    set unit_pivot [cheap_vecnorm $pivot_vec]
    set rho [vecsub $r [vecscale [vecdot $r $unit_pivot] $unit_pivot]] ;# find the vector rejection in order to get the "arm" length
    set rho_sq [ vecdot $rho $rho ]
    set moment [ expr "$moment + ( $mass * $rho_sq )" ]
  }
  return $moment
}

proc get_inertia_matrix_slow { weights coords COM } { ;# will need to replace this function with an equivalent one in C
  set n [llength $weights]
  #set inertia_matrix { { 0.0 0.0 0.0 } { 0.0 0.0 0.0 } { 0.0 0.0 0.0 } }
  set Ixx 0.0; set Iyy 0.0; set Izz 0.0
  set Ixy 0.0; set Ixz 0.0; set Iyz 0.0
  for {set i 0} {$i < $n} {incr i} {
    set m [lindex $weights $i]
    set x [expr "[lindex [lindex $coords $i] 0] - [lindex $COM 0]"]
    set y [expr "[lindex [lindex $coords $i] 1] - [lindex $COM 1]"]
    set z [expr "[lindex [lindex $coords $i] 2] - [lindex $COM 2]"]
    set Ixx [expr "$Ixx + $m * ($y * $y + $z * $z)"]
    set Iyy [expr "$Iyy + $m * ($x * $x + $z * $z)"]
    set Izz [expr "$Izz + $m * ($x * $x + $y * $y)"]
    set Ixy [expr "$Ixy - $m * ($x * $y)"]
    set Ixz [expr "$Ixz - $m * ($x * $z)"]
    set Iyz [expr "$Iyz - $m * ($y * $z)"]
  }
  set inertia_matrix " {$Ixx $Ixy $Ixz} {$Ixy $Iyy $Iyz} {$Ixz $Iyz $Izz} "
  return $inertia_matrix
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
set CROSS_QUAT $MAIN_NEIGHBOR
set anchor_to_cross [quat_mult $CROSS_QUAT [quat_inv $ANCHOR_QUAT]]
set anchor_to_anticross [quat_inv $anchor_to_cross]
set ANTICROSS_QUAT [ quat_to_anchor_hemisphere  $ANCHOR_QUAT [quat_mult $anchor_to_anticross $ANCHOR_QUAT]]
output "CROSS_QUAT: $CROSS_QUAT"
output "ANTICROSS_QUAT: $ANTICROSS_QUAT"
# construct the sampling plane
set EQ_ANGLE [get_quat_angle $MAIN_NEIGHBOR $ANCHOR_QUAT]
set NEIGHBOR_QUATS [quats_to_anchor_hemisphere $ANCHOR_QUAT $NEIGHBOR_QUATS] ;# fix the neighbors so they are all in the same hemisphere as the anchor

set anchor_to_main_neighbor_quat [quat_mult $MAIN_NEIGHBOR [quat_inv $ANCHOR_QUAT]] ;# the normal points out from the plane
output "main_neighbor: $MAIN_NEIGHBOR"
output "anchor_to_main_neighbor_quat: $anchor_to_main_neighbor_quat"
# construct the list that holds all the quaternions that point from the neighbors to the anchor
#set neighbor_to_anchor_quats ""
#set neighbor_to_anchor_angles ""
set NUM_NEIGHBORS [llength $NEIGHBOR_QUATS]
set max_quat_dot 0.0
foreach neighbor_quat $NEIGHBOR_QUATS {
  if { [quat_dot $neighbor_quat $ANCHOR_QUAT] >= 1.0 } {continue} ;# then this is ourself, so skip
  #lappend neighbor_to_anchor_quats [quat_mult [quat_inv $ANCHOR_QUAT] $neighbor_quat]
  #lappend neighbor_to_anchor_angles [get_quat_angle $neighbor_quat $ANCHOR_QUAT ]

  if { [quat_dot $neighbor_quat $ANCHOR_QUAT] > $max_quat_dot} {
    set max_quat_dot [vecdot $neighbor_quat $ANCHOR_QUAT]
    set closest_quat_angle [expr "acos($max_quat_dot)"] ;# not multiplied by 2.0 because it's halfway between
  }
}
#output "neighbor_to_anchor_angles: $neighbor_to_anchor_angles"
output "NUM_NEIGHBORS: $NUM_NEIGHBORS"
output "max_quat_dot: $max_quat_dot"
output "closest_quat_angle: $closest_quat_angle"

# calculate the quaternion potential coefficients based on the CROSS and the ANTICROSS

set QUAT_POTENTIAL_INDECES {0 1 2 3} ;# the indeces of the ligand quaternion that determine the potential energy
set QUAT_POTENTIAL_COEFFS [cheap_vecnorm [vecsub $CROSS_QUAT $ANTICROSS_QUAT] ] ;# the coefficients are simply the difference between the CROSS and the ANTICROSS

output "QUAT_POTENTIAL_COEFFS: $QUAT_POTENTIAL_COEFFS"

proc calcforces { } {
  global ligatoms; global recatoms; global ANCHOR_QUAT; global CROSS_QUAT; global ANTICROSS_QUAT; global FORCE_CONSTANT; global NULL_QUAT;
  global EQ_ANGLE; global ligrefcoords; global ligweights; global recrefcoords; global recweights; global closest_quat_angle;
  global NEIGHBOR_QUATS; global RECROT; global ligrefquat; global ligquat_addition; global ligrefCOM; global recrefCOM; global recrefquat
  global old_ligquat; global neighbor_to_anchor_quats; global neighbor_to_anchor_angles; global old_d_theta;
  global QUAT_POTENTIAL_INDECES; global QUAT_POTENTIAL_COEFFS
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
    set ligrefquat_raw [confs_to_quat $ligrefcoords $ligcoords $ligrefCOM $ligrefCOM $ligweights] ;# should be NULL_QUAT
    set ligquat_addition [quat_mult $ANCHOR_QUAT [quat_inv $ligrefquat_raw]] ;# the quaternion that points from the reference ligand quaternion to ANCHOR_QUAT
    if {$RECROT == "True"} {
      set recrefquat [confs_to_quat $recrefcoords $reccoords $recrefCOM $recrefCOM $recweights] ;# should be NULL_QUAT
    } else {
      set recrefquat $NULL_QUAT
    }
    set ligrefquat [quat_to_anchor_hemisphere $ANCHOR_QUAT [quat_mult $ligquat_addition [quat_mult $ligrefquat_raw [quat_inv $recrefquat]]]] ;# factor out the effect that the receptor has rotated as well as ANCHOR_QUAT
    set EQ_ANGLE [get_quat_angle $ligrefquat $CROSS_QUAT]
    set old_ligquat $ligrefquat
    set old_d_theta 0.0
    output "recrefquat: $recrefquat"
    output "ligrefquat: $ligrefquat"
    output "ligquat_addition: $ligquat_addition"


  }

  #set ligCOM [getCOM_fast $ligweights $ligcoords ] ;# WARNING: getCOM_fast has a bug
  #set recCOM [getCOM_fast $recweights $reccoords ] ;# WARNING: getCOM_fast has a bug
  set ligCOM [getCOM $ligcoords $ligweights ]
  set recCOM [getCOM $reccoords $recweights ]

  #output "ligCOM: $ligCOM"
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
  output "ligquat: $ligquat"
  if { [vecdot $ligquat $old_ligquat] < 0.9} {
    output "LARGE QUATERNION JUMP! [vecdot $ligquat $old_ligquat]"
    output "ligquat: $ligquat"
    output "old_ligquat: $old_ligquat"
    #abort
  }

  # impose 1st force: maintaining a desired dot between cross_quat and anticross_quat
  #set lig_torque1 [get_torque_from_quat_potential $ligquat $FORCE_CONSTANT $QUAT_POTENTIAL_INDECES $QUAT_POTENTIAL_COEFFS ] ;# this form doesn't seem to be working properly
  set q_prime1 [quat_mult [quat_inv $CROSS_QUAT] $ligquat ]
  set lig_torque1 [get_torque_from_quat_potential_angle $q_prime1 $EQ_ANGLE $FORCE_CONSTANT ]
  #output "lig_torque1: $lig_torque1"
  #set crossdot [vecdot $ligquat $CROSS_QUAT] ;# get dot product of two quaternions
  #set anticrossdot [vecdot $ligquat $ANTICROSS_QUAT]
  #output "crossdot: $crossdot"
  #output "anticrossdot: $anticrossdot"
  #output "difference: [expr {$crossdot - $anticrossdot}]"
  set anchor_quat_angle [get_quat_angle $ligquat $ANCHOR_QUAT]
  set main_neighbor_quat_angle [get_quat_angle $ligquat $CROSS_QUAT]
  # impose 2nd force: maintain a minimum quat distance from ANCHOR_QUAT
  set lig_torque2 {0 0 0}
  #output "cross quat angle: $cross_quat_angle"
  output "anchor_quat_angle $anchor_quat_angle"
  output "main_neighbor_quat_angle: $main_neighbor_quat_angle"
  output "EQ_ANGLE: $EQ_ANGLE"
  output "closest_quat_angle $closest_quat_angle"
  if { [expr "abs($anchor_quat_angle)"] > $closest_quat_angle } {
    set q_prime2 [quat_mult [quat_inv $ANCHOR_QUAT] $ligquat ]
    set lig_torque2 [get_torque_from_quat_potential_angle $q_prime2 $closest_quat_angle $FORCE_CONSTANT]
    output "lig_torque2: $lig_torque2"
  }
  #if { ( [expr "abs([lindex $ligquat 0])"] < 0.01 ) } {
  #    output "LIGQUAT COLLINEAR! $ligquat"
  #    #abort
  #  }
  # add the torques
  set lig_torque [vecadd $lig_torque1 $lig_torque2 ]
  # convert torque to a linear force; impose that linear force on all atoms in the ligand
  if {[veclength $lig_torque] > 0.000001 } { ;# if the torque vector is not tiny
    set ang_accel [get_angular_accel $lig_torque $ligweights $ligCOM $ligcoords]
    foreach atom $ligatoms {
      set r [vecsub $coords($atom) $ligCOM]
      set unit_ang_accel [cheap_vecnorm $ang_accel]
      set rho [vecsub $r [vecscale [vecdot $r $unit_ang_accel] $unit_ang_accel]] ;# find the vector rejection in order to get the "arm" length
      set force [ vecscale [ veccross $ang_accel $rho ] $masses($atom)]
      addforce $atom $force
    }
  } else {
    set force {0 0 0}
  }
  set old_ligquat $ligquat
  return
}
