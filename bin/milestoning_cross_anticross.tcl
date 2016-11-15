# Variables/Constants required from NAMD input file
#  ID ;# the ID of this trajectory
#  SCRIPT_INTERVAL ;# number of timesteps before script should be evaluated. Examples: 1=script evaluated every timestep. 5=script eval'd every 5th timestep, ...
#  LIGROT ;# whether we care about ligand rotation
#  ABORT_ON_CROSSING
#  PHASE  ;# can be 'forward' or 'reverse'
#  CARE_ABOUT_SELF ;# whether we care about the milestone's self in this simulation
#  MILESTONE_LIST ;# a list of milestones to be converted to arrays
#  whoami  ;# the index of the current
#  GRID_EDGE_RAD  ;# radius to the edge of the grid from the center
#  LIGRANGE ;# the indeces of the atoms that define the ligand
#  RECRANGE ;# the indeces of the atoms that define the receptor
#  RECPA1_LIST ;# the list of relatively immobile atoms that define the receptor's 1st PA
#  RECPA3_LIST ;# the list of relatively immobile atoms that define the receptor's 3rd PA
#  FHPD_FILE ;# the location of the file where we will write FHPD trajectory ID's


set outputprefix "SEEKR: "
set first_timestep_tolerance 10 ;# the number of timesteps to wait before aborting the trajectory during the reversal phase
set first_milestone_hit False ;# whether any milestone has been hit since the start of the simulation
set MAX_NUM_STEPS 100000

if {!([info exists ID])} {set ID 0}
if {!([info exists SITEID])} {set SITEID 0}
#set LA_src "$TEMPLATE_LA_src" ;#NOTE: this will need to be defined by the user
#if {$LA_src != ""} { source $LA_src } ;# import the tcl linear algebra package

set incubation_start 0

set NULL_QUAT "1 0 0 0"
set PI 3.1415926

proc add_id_to_FHPD {id} {
  ###########################################################################################
  # Opens a file named $FHPD_FILE and writes integer $id to it
  ###########################################################################################
  global FHPD_FILE
  # opens a file to write the current simulation's id to be used later in the forward phase
  if { [file exists $FHPD_FILE] } {
    set ourfile [open $FHPD_FILE a] ;# then just append to it
  } else {
    set ourfile [open $FHPD_FILE w] ;# then just append to it
  }
  puts $ourfile $id
  close $ourfile
}

proc switch_grid_point {new_id} {
  # given a new index for the planes list, will construct a new list of neighbors
  global whoami
  set whoami $new_id
}

proc switch_rot_point {new_id} {
  # given a new index for the planes list, will construct a new list of neighbors
  global whoami_rot
  set whoami_rot $new_id
}

#proc closest_grid_point {ligcenter neighbors} {
#  ###########################################################################################
#  # finds point in $neighbors closest to $ligcenter.
#  ###########################################################################################
#  array set array_neighbors $neighbors ;# pass argument as a list, convert to array once we're here
#  set closestdist 999999
#  set closestneighbor None
#  foreach neighbor [array names array_neighbors] { ;# loop through each neighbor to see which one we are closest to
#    #puts "neighbor: $neighbor, $array_neighbors($neighbor)"
#    set neighbor_lig_dist [cheap_vecdist $neighbor $ligcenter]
#    if {$neighbor_lig_dist < $closestdist} { ;# set this to be the closest neighbor
#      set closestdist $neighbor_lig_dist
#      set closestneighbor $array_neighbors($neighbor)
#    }
#  }
#  return $closestneighbor
#}

proc output { message } { ;# automatically appends caption to stdout from script
  global outputprefix;
  puts "$outputprefix$message"
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


###########################################################################################
# MATH: missing vector/matrix math functions
###########################################################################################
proc cheap_vecdist {vec1 vec2} { ;# distance between the ends of two vectors
  return [veclength [vecsub $vec2 $vec1]]
}

proc cheap_vecnorm {a} { ;# normalize a vector to be of length 1.0
  set scalefactor [expr "1/[veclength $a]"]
  return [vecscale $a $scalefactor]
}

###########################################################################################
# INITIALIZATIONS: continued
###########################################################################################

set liglist [parse_selection $LIGRANGE]
set reclist [parse_selection $RECRANGE]

set ligatoms {}
foreach atomid $liglist {
  lappend ligatoms $atomid
}
set recatoms {}
foreach atomid $reclist {
  lappend recatoms $atomid
}
foreach atom $ligatoms {
  addatom $atom
}
foreach atom $recatoms {
  addatom $atom
}
set MILESTONE_LIST_COPY $MILESTONE_LIST
set MILESTONE_LIST ""
foreach milestone_unset $MILESTONE_LIST_COPY {
  array set milestone $milestone_unset
  if {(($milestone(shape) == "sphere" ) || ($milestone(shape) == "plane" )) && ($milestone(center_type) == "atom")} { ;# then we will want to add these atoms to a coordinate group
    set group [addgroup $milestone(center_indeces)] ;# define the group that represents the center of this milestone
    #foreach center_index $milestone(center_indeces) {
    #  addatom $center_index
    #}
    set milestone(center_com_id) $group ;# NOTE: we need to base this off a true center of mass, not atomic group

  } elseif { ($milestone(shape) == "rotational") && ($LIGROT == "True") } {
    if {$milestone(id) == $whoami_rot} {
      set ANCHOR_QUAT $milestone(quaternion)
    }
  }

  set milestone_unset [array get milestone]
  lappend MILESTONE_LIST $milestone_unset
}

#set recidlist {}
#foreach atom $reclist {
#  lappend recidlist [addgroup $atom]
#}
#set PA1idlist {}
#foreach atom $RECPA1_LIST {
#  lappend PA1idlist [addgroup $atom]
#}
#set recPA1id [addgroup $RECPA1_LIST]
#set PA3idlist {}
#foreach atom $RECPA3_LIST {
#  lappend PA3idlist [addgroup $atom]
#}
#set recPA3id [addgroup $RECPA3_LIST]

#set PA1ligidlist {}
#foreach atom $LIGPA1_LIST {
#  lappend PA1ligidlist [addgroup $atom]
#}
#set ligPA1id [addgroup $LIGPA1_LIST]
#set PA3ligidlist {}
#foreach atom $LIGPA3_LIST {
#  lappend PA3ligidlist [addgroup $atom]
#}
#set ligPA3id [addgroup $LIGPA3_LIST]

#set startrecPA {}
#set startligPA {}
#set gridcenter [vecsub $ligCOM $recCOM] ;# vector from rec COM with rec PAs as unit vector describing the center of the grid


array set test_milestone1 {
  id 0
  shape sphere
  center_type atom
  absolute False
  center_indeces { 613 }
  center_com_id ""
  radius 10.0
  side inside
  end False
} ;# NOTE: center_com_id is supposed to start empty; we fill it later

array set test_milestone2 {
  id 1
  shape plane
  center_type coord
  absolute True
  center {12.3 45.6 78.9}
  normal {0.0 0.0 1.0}
  side negative
}



proc check_spherical_milestone {coord1 coord2 radius } {
  ###########################################################################################
  # function returns 'inside' if trajectory on inside of sphere, 'outside' otherwise
  ###########################################################################################
  set rad_sq [ expr "$radius * $radius" ]
  if {[ llength $coord1 ] != 3 } {error "check_spherical_milestone: invalid coord1 value (must be list of length 3): $coord1"}
  if {[ llength $coord2 ] != 3 } {error "check_spherical_milestone: invalid coord2 value (must be list of length 3): $coord2"}
  set diff [vecsub $coord2 $coord1]
  set dist_sq [vecdot $diff $diff]
  set dist [expr "sqrt($dist_sq)"]
  #print "dist: $dist radius: $radius"
  if { $dist_sq < $rad_sq } { return inside } else { return outside }
}

proc check_planar_milestone {coord milestone_center milestone_normal } {
  ###########################################################################################
  # function returns 'negative' if trajectory on negative side of plane relative to the normal, 'positive' otherwise
  ###########################################################################################
  set coord_center_diff [vecsub $coord $milestone_center] ;# the vector that points from the milestone center to the trajectory coord
  set coord_normal_dot [vecdot $coord_center_diff $milestone_normal] ;# the dot between the diff and normal
  if {$coord_normal_dot > 0} { ;# then the trajectory and the normal are on the same side of the milestone
    return positive
  } else {
    return negative
  }
}

proc make_reversals { flag } {
  # This function checks if velocity reversals have been properly implemented in this NAMD.
  # If so, it flags NAMD to reverse velocities in the current step
  output "REVERSING ALL VELOCITIES!!!"
  if {[ catch {setboundaryflag $flag} ] != 0} {output "ERROR: setboundaryflag command not recognized in this compilation of NAMD."; abort}
}

proc get_quat_angle { quat1 quat2 } { ;# find the angle between two quaternions
  # WARNING: does not use quat_dot, so the quaternions MUST be in the proper hemisphere
  set  [vecdot $quat1 $quat2] ;# get dot product of two quaternions
  return [expr "2.0 * acos($quatdot)"] ;# arccos of the dot product is the angle
}

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
  return $result
}

proc quats_to_torque {quat1 quat2 eq_theta force_coeff} {
  # given two quaternions, will calculate the torque needed to maintain a constant equilibrium angle for quat1 to reach quat2
  set q_prime [quat_mult $quat2 [quat_inv $quat1]] ;# this is the quaternion to make quat1 rotate to quat2
  set axis [lrange $q_prime 1 end] ;# the axis of rotation between the two
  if {[veclength $axis] <= 0.00000001} {return {0.0 0.0 0.0}}
  set theta [expr "2.0 * acos([lindex $q_prime 0])"] ;# the angle between the two quats
  set d_theta [expr "$theta - $eq_theta"] ;# the difference between the quat angle and the equilibrium angle
  set torque_magnitude [expr "-2.0 * $force_coeff * $d_theta"] ;# magnitude of the torque
  set torque_unit [cheap_vecnorm $axis] ;# the direction of the torque
  return [ vecscale $torque_unit $torque_magnitude] ;# scale the direction by the magnitude to get the torque vector
}

proc quat_to_torque {q_prime theta eq_theta force_coeff} {
  # given one quaternion, will calculate the torque needed to maintain a constant equilibrium angle for quat1 to reach quat2
  #set q_prime [quat_mult $quat2 [quat_inv $quat1]] ;# this is the quaternion to make quat1 rotate to quat2
  set axis [lrange $q_prime 1 end] ;# the axis of rotation between the two
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

proc calculate_anticross { refquat crossquat  } {
  # given a reference quaternion and a cross quaternion, will return an anticross quaternion
  set anchor_to_cross [quat_mult $crossquat [quat_inv $refquat]]
  set anchor_to_anticross [quat_inv $anchor_to_cross]
  return [ quat_to_anchor_hemisphere $refquat [quat_mult $anchor_to_anticross $refquat]] ;# return the anticross
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

proc rotate_by_quat { rot_quat start_point { newCOM "0 0 0" } { oldCOM "0 0 0" } } {
  # given a quaternion, a start_point, and optionally, an old Center-Of-Mass and a new Center-Of-Mass,
  #  will translate the point from the old COM to the new COM, and then rotate by the quaternion
  set orig_point [vecsub $start_point $oldCOM] ;# first move the point to the origin
  set angle [expr "2.0 * acos([lindex $rot_quat 0])"] ;# the angle of rotation about the axis in radians
  set axis [lrange $rot_quat 1 end] ;# the axis of rotation
  if { $angle == 0.0 } {return $start_point} ;# if there's no rotation then return the original point
  set rotmat [transabout $axis $angle rad] ;# find the rotation matrix about the quaternion
  set new_orig_point [vectrans $rotmat $orig_point] ;# rotate point about the origin
  return [vecadd $new_orig_point $newCOM] ;# translate to the new coordinates and return
}

#puts "testing transabout(rad): [transabout {0 1 0} 1.0 rad]"
#puts "testing transabout(deg): [transabout {0 1 0} 114.592 deg]"

proc calcforces { } {

  global neighbors; global ligatoms; global recatoms; global ligweights; global recweights; global ligrefcoords; global recrefcoords
  global raw_neighbors; global gridcenter; global whoami; global whoami_rot; global ligquat_addition; global ligrefCOM; global recrefCOM
  global incubation_start; global recCOM; global SCRIPT_INTERVAL; global PHASE; global CARE_ABOUT_SELF; global LIGROT; global RECROT
  global incubation_start_rot; global first_milestone_hit_rot; global oldquat; global newquat
  global ABORT_ON_CROSSING; global LIGRANGE; global RECRANGE; global ANCHOR_QUAT; global NULL_QUAT
  global FHPD_FILE; global MILESTONE_LIST; global first_timestep_tolerance; global ID; global first_milestone_hit;
  global SITEID; global ligrefquat; global recrefcoords; global ligquat; global MAX_NUM_STEPS

  set stepnum [getstep] ;# step number
  setboundaryflag 0 ;# make sure the reversals are turned on only for the steps we want them turned on for
  if { [expr "$stepnum % $SCRIPT_INTERVAL"] != 0 } {return} ;# skip evaluation if not on an evaluation timestep

  # get coordinates of important atom groups
  loadcoords coords
  loadmasses masses
  set ligcoords [getcoords $ligatoms]
  set reccoords [getcoords $recatoms]

  # everything to perform on the first step
  if { $stepnum == 0 } {
    if {$PHASE == "reverse"} { make_reversals 1 } ;# for the reversal phase of the milestoning: reverse all velocities to start with. NOTE: requires special build of NAMD
    switch_grid_point $whoami
    if {$LIGROT == "True"} {
      set ligrefcoords $ligcoords
      set recrefcoords $reccoords
      set ligweights [getweights $ligatoms]
      set recweights [getweights $recatoms]
      set ligrefCOM [getCOM $ligcoords $ligweights]
      set recrefCOM [getCOM $reccoords $recweights]
      set ligrefquat_raw [confs_to_quat $ligrefcoords $ligrefcoords $ligrefCOM $ligrefCOM $ligweights] ;# should be NULL_QUAT
      set ligquat_addition [quat_mult $ANCHOR_QUAT [quat_inv $ligrefquat_raw]] ;# the quaternion that points from the reference ligand quaternion to ANCHOR_QUAT
      set recrefquat [confs_to_quat $recrefcoords $reccoords $recrefCOM $recrefCOM $recweights] ;# should be NULL_QUAT
      set ligrefquat [quat_to_anchor_hemisphere $ANCHOR_QUAT [quat_mult $ligquat_addition [quat_mult $ligrefquat_raw [quat_inv $recrefquat]]]] ;# factor out the effect that the receptor has rotated as well as ANCHOR_QUAT
      set newquat $ligrefquat; set oldquat $ligrefquat
    }
    output "whoami: $whoami"
    output "whoami_rot: $whoami_rot"
    output "MILESTONES"
    set counter 0
    foreach i $MILESTONE_LIST {
      output "milestone: $i"
      array set milestone $i
      #output "recrefCOM: $recrefCOM"
      if { [info exists milestone(normal)] == 1 } { set milestone(center) [ vecadd $recrefCOM [ vecscale $milestone(normal) $milestone(distance) ] ] }
      set milestone_unset [array get milestone]
      #output "milestone_unset: $milestone_unset"
      set MILESTONE_LIST [lreplace $MILESTONE_LIST $counter $counter $milestone_unset]
      incr counter
    }
    #for {set i 0} {$i < [expr "[llength $raw_neighbors]/2"]} {incr i} {lappend neighbors [vecadd [lindex $raw_neighbors [expr "$i*2"]] $gridcenter]; lappend neighbors [lindex $raw_neighbors [expr "$i*2 + 1"]]}
  } ;# vector from rec COM with rec PAs as unit vector describing the center of the grid

  # EVERY TIMESTEP (AT INTERVAL) EXECUTES CODE BELOW

  set counter 0
  set crossed "False"; set crossed_rot "False"; set aborting False

  if {$LIGROT == "True"} { ;# find general information that will apply to any milestone
    set ligCOM [getCOM_fast $ligweights $ligcoords ]
    set recCOM [getCOM_fast $recweights $reccoords ]
    set ligquat_raw [confs_to_quat $ligcoords $ligrefcoords $ligCOM $ligrefCOM $ligweights] ;# use Quaternion Least Squares to obtain the rotation quaternion
    if {$RECROT == "True"} { ;# if there is a receptor specified and its rotations are turned 'on'
      set recquat [quat_to_anchor_hemisphere $ANCHOR_QUAT [confs_to_quat $reccoords $recrefcoords $recCOM $recrefCOM $recweights]] ;# get the rotation that recquat has made
    } else { ;# then we are dealing with absolute ligand rotation: it is not rotating with the receptor
      set recquat $NULL_QUAT
    }
    set ligquat [quat_to_anchor_hemisphere $ANCHOR_QUAT [quat_mult $ligquat_addition [quat_mult $ligquat_raw [quat_inv $recquat]]]] ;# factor out the effect that the receptor has rotated as well as ANCHOR_QUAT
    output "ligquat: $ligquat"
  }
  # handle each milestone
  foreach milestone_unset $MILESTONE_LIST {
    array set milestone $milestone_unset ;# first we have to unpack the milestone array
    #if { (($CARE_ABOUT_SELF == False) || ($stepnum <= $first_timestep_tolerance)) && ($whoami == $milestone(id))} { incr counter; continue } ;# if milestones are ignoring self-crossings, then skip this milestone. Also skip if the milestone is crossing itself, and the trajectory needed some time to move to a side. Increment the counter anyways

    # check milestone center_type (how the milestone moves with the receptor)
    if {(($milestone(shape) == "sphere" ) || ($milestone(shape) == "plane" ))} {
      if { $milestone(center_type) == "atom" } { ;# then the milestone is centered on an atom group
        set milestone_center_id $milestone(center_com_id)
        set milestone_center $coords($milestone_center_id)
        if {[ info exists milestone(normal) ]} {set milestone_normal $milestone(normal)} ;# only valid for planar milestones
      } elseif { $milestone(center_type) == "coord" } { ;# then the coordinate is relative to the rec COM and PA's
        # find milestone center/normal
        #output "coord type milestones not yet implemented"
        if { $milestone(absolute) == True } { ;# then the milestone doesn't move with the receptor
          set milestone_center $milestone(center)
          if {[ info exists milestone(normal) ]} {set milestone_normal $milestone(normal)} ;# only valid for planar milestones
        } else { ;# then the milestones DO move with the receptor
          # calculate milestone coordinates relative to the mobile receptor
          if { $RECROT == "True" } {
            set recquat [confs_to_quat $recrefcoords $reccoords $recCOM $recrefCOM $recweights] ;# get the rotation that recquat has made
          } else {
            set recquat $NULL_QUAT
          }
          # move the milestone center relative to the rotations
          set milestone_center [rotate_by_quat $recquat $milestone(center) $recCOM $recrefCOM]
          #set milestone_center [ calc_ref_movements $milestone(center) $recCOM $recrefCOM $recPA $startrecPA ]
          #output "recCOM: $recCOM, recPA: $recPA, milestone_center: $milestone_center, milestone(center): $milestone(center)"
          if {[ info exists milestone(normal) ]} { ;# with planar milestones
            #set milestone_normal_plus_center [ calc_ref_movements [vecadd $milestone(center) $milestone(normal)] $recCOM $recrefCOM $recPA $startrecPA ] ;# in order to get the normal, need to add normal to center, calculate shift in frame of ref.
            #set milestone_normal [vecsub $milestone_normal_plus_center $milestone_center] ;# ... and then subtract away the new center to get the new normal
            # set the new milestone normal so that we can keep track of how it has changed over time
            set milestone_normal_plus_center [rotate_by_quat $recquat [vecadd $milestone(center) $milestone(normal) ] $recCOM $recrefCOM]
            set milestone_normal [vecsub $milestone_normal_plus_center $milestone_center] ;# ... and then subtract away the new center to get the new normal
          }
        }
      }
    }
    # check the milestone shape
    if { $milestone(shape) == "sphere" } {
      set result [check_spherical_milestone $ligCOM $milestone_center $milestone(radius)]
    } elseif { $milestone(shape) == "plane" } {
      set result [check_planar_milestone $ligCOM $milestone_center $milestone_normal]
    } elseif {( $milestone(shape) == "rotational") && ($LIGROT == True)} {
      # check rotational milestones for each milestone
      output "WARNING: LIGAND ROTATIONS NOT COMPLETELY IMPLEMENTED!!"
      # find whether ligand has crossed this current milestones
      # if the ligand is close to the this milestone's cross then the side is "cross", if closer to the anticross then the side is "anticross"
      # may need to differentiate between whether this is the source milestone being monitored for the reversal phase
      #  or an adjacent milestones being monitored for a different type of crossing event
      if {($whoami_rot == $milestone(id)) && ($CARE_ABOUT_SELF == True)} { ;# then we are monitoring our current milestone for crossings
        # first define all cross/anticross pairs
        # then define the side of each surface the ligand is on
        foreach cross_anticross_pair $milestone(cross_anticross_list) {
          set cross [lindex $cross_anticross_pair 0]
          set anticross [lindex $cross_anticross_pair 1]
          if {[quat_dot $cross $ligquat] > [quat_dot $anticross $ligquat]} {
            set result "cross"
          } else {
            set result "anticross"
          }

        }
      } else {
        # then figure out which cross/anticross pair we are using to define the crossing surface...
        set cross $ANCHOR_QUAT
        set anticross [calculate_anticross $milestone(quaternion) $cross] ;# NOTE: rewrite so this doesn't have to be done each timestep
        if {[quat_dot $cross $ligquat] > [quat_dot $anticross $ligquat]} {
          set result "cross"
        } else {
          set result "anticross"
        }
      }

    } ;# elseif other milestone types ???

    # check for a crossing event
    if {[info exists milestone(siteid)]} { ;# if in this system, the milestones contain siteids
      set milestone_siteid $milestone(siteid)
    } else { # then it hasn't been implemented yet
      set milestone_siteid 0
    }
    if { ($result != $milestone(side)) && ($milestone_siteid == $SITEID) } { ;# then a crossing event has occurred
      output "crossing event detected for milestone: $milestone(id), result: $result"
      output "Milestones: $MILESTONE_LIST"
      set milestone(side) $result ;# set the new side that the trajectory is on
      if {$stepnum > $first_timestep_tolerance && ($whoami != $milestone(id) || $CARE_ABOUT_SELF == True)} { # then this is a nonzero step and the milestone sides have already been set
        if { $milestone(shape) == "rotational" } {
          set crossed_rot $milestone(id) ;# say that we have crossed a milestone
          set oldquat $newquat
          set newquat $milestone(quaternion)
        } else { ;# then its not a rotational milestone
          set crossed $milestone(id) ;# say that we have crossed a milestone
        }
        if { ($PHASE == "forward") && ($first_milestone_hit == "False") && ($stepnum > 0) } {
          if { (($milestone(shape) != "rotational") && ($whoami != $milestone(id))) || (($milestone(shape) == "rotational") && ($whoami_rot != $milestone(id))) } {
            # then we have a problem, because the forward phase is supposed to hit the $whoami milestone first before anything else
            output "FORWARD ERROR: Forward phase trajectory crossed different milestone that the starting milestone. Wanted: $milestone(id), crossed: $whoami. This trajectory must be discarded. Aborting..."
            abort
          } else {
            # then the forward phase has crossed the starting plane, as expected.
            output "Forward phase crossed first positional milestone. Resetting incubation time which was at: [expr $stepnum - $incubation_start]"
            set first_milestone_hit True
            set incubation_start $stepnum
           }
        }

        if { $milestone(end) == True } { set aborting True } ;# if this is an ending milestone, then abort simulation
      }
      output "Milestones: $MILESTONE_LIST"
    }
    # replace the member of the list with the milestone array
    set milestone_unset [array get milestone]
    set MILESTONE_LIST [lreplace $MILESTONE_LIST $counter $counter $milestone_unset]
    incr counter
  }

  if {$crossed != "False"} { ;# print useful info and, if necessary, abort simulation
    output "Milestones: $MILESTONE_LIST"
    output "GRID TRANSITION: source: $whoami, destination: $crossed, current step: $stepnum, ligand COM: $ligCOM, receptor COM: $recCOM, receptor start COM: $recrefCOM, incubation time: [expr {$stepnum - $incubation_start}] timesteps"
    if { $PHASE == "reverse" && $whoami != $crossed } { output "adding to FHPD file"; add_id_to_FHPD $ID } ;# then the reversal has crossed another milestone, and we should record it
    set incubation_start $stepnum
    switch_grid_point "$crossed"
    if { ($ABORT_ON_CROSSING == True) || ($aborting == True) || (($stepnum > $MAX_NUM_STEPS) && ($MAX_NUM_STEPS != 0))} { abort }
  }
  if {$crossed_rot != "False"} { ;# print useful info and, if necessary, abort simulation

      output "ROTATION TRANSITION: source: $whoami_rot, destination: $crossed_rot, current step: $stepnum, ligand COM: $ligCOM, receptor COM: $recCOM, receptor start COM: $recrefCOM, incubation time: [expr {$stepnum - $incubation_start}] timesteps, ligquat: $ligquat, oldquat: $oldquat, newquat: $newquat"
      if { $PHASE == "reverse" && $whoami_rot != $crossed_rot } { output "adding to FHPD file"; add_id_to_FHPD $ID } ;# then the reversal has crossed another milestone, and we should record it
      set incubation_start $stepnum
      switch_rot_point "$crossed_rot"
      if { ($ABORT_ON_CROSSING == True) || ($aborting == True) || (($stepnum > $MAX_NUM_STEPS) && ($MAX_NUM_STEPS != 0))} { abort }
    }
}

###########################################################################################
# optimization section, variables that have been set to optimize simulation
###########################################################################################


