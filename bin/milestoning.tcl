# Variables/Constants required from NAMD input file (NOTE: need to verify these...)
#  ID ;# the ID of this trajectory
#  VEL_ID ;# the velocity ID of this trajectory, in case the same position conformation was relaunched many times with new velocities.
#  SCRIPT_INTERVAL ;# number of timesteps before script should be evaluated. Examples: 1=script evaluated every timestep. 5=script eval'd every 5th timestep, ...
#  ABORT_ON_CROSSING
#  STAGE  ;# can be 'forward' or 'reverse'
#  CARE_ABOUT_SELF ;# whether we care about the milestone's self in this simulation
#  MILESTONE_LIST ;# a list of milestones to be converted to arrays
#  whoami  ;# the index of the current milestone state
#  GRID_EDGE_RAD  ;# radius to the edge of the grid from the center
#  LIGRANGE ;# the indeces of the atoms that define the ligand
#  RECRANGE ;# the indeces of the atoms that define the receptor
#  PERFECT_STARTING_DIST ;# whether the system has been truly constrained to the milestone in the eq. sims. Should be False for umbrella sampling

set outputprefix "SEEKR: "
set first_timestep_tolerance 10 ;# the number of timesteps to wait before aborting the trajectory during the reversal stage
set first_milestone_hit False ;# whether any milestone has been hit since the start of the simulation
set MAX_NUM_STEPS 1000000 ;# maximum number of steps before aborting
#set incubation_start 0
if {!([info exists ID])} {set ID 0}
if {!([info exists VEL_ID])} {set VEL_ID 0}
if {!([info exists SITEID])} {set SITEID 0}
if {!([info exists PERFECT_STARTING_DIST])} {set PERFECT_STARTING_DIST "False"} ;# set it False by default
set PI 3.1415926

proc switch_grid_point {new_id} {
  # given a new index for the planes list, will construct a new list of neighbors
  global whoami
  set whoami $new_id
}

proc outprint { message } { ;# automatically appends caption to stdout from script
  global outputprefix;
  puts "$outputprefix$message"
}

proc range { start cutoff finish {step 1} } {
  # gives a range of numbers from start to finish, incrementing by step
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

  }
  set milestone_unset [array get milestone]; array unset milestone
  lappend MILESTONE_LIST $milestone_unset
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
  #set dist [expr "sqrt($dist_sq)"]
  #print "dist: $dist"
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

# start all the ligand information
proc getCOM { atoms weights } { ;# can this be removed?
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

proc append_to_file { filename content } { ;# quick write to the end of a file
  set ourfile [open $filename a]
  puts $ourfile $content
  close $ourfile
}

proc calcforces { } {
  global neighbors; global ligatoms; global recatoms; global ligweights; global recweights; global ligrefcoords; global recrefcoords
  global raw_neighbors; global gridcenter; global whoami; global whoami_rot; global ligquat_addition; global ligrefCOM; global recrefCOM
  global incubation_start; global recCOM; global SCRIPT_INTERVAL; global STAGE; global CARE_ABOUT_SELF; global LIGROT; global RECROT
  global incubation_start_rot; global first_milestone_hit_rot; global oldquat; global newquat; global EQ_ANGLE
  global ABORT_ON_CROSSING; global LIGRANGE; global RECRANGE; global ANCHOR_QUAT; global NULL_QUAT
  global MILESTONE_LIST; global first_timestep_tolerance; global ID; global VEL_ID; global first_milestone_hit;
  global SITEID; global ligrefquat; global recrefcoords; global ligquat; global MAX_NUM_STEPS
  global PERFECT_STARTING_DIST; global reverse_failed; global aborting; global stepnum; global FWD_FILENAME; global REV_FILENAME; global RESTART_FREQ
  #set stepnum [getstep] ;# step number but this isn't right
  #setboundaryflag 0 ;# make sure the reversals are turned on only for the steps we want them turned on for
  if { [expr "[getstep] % $SCRIPT_INTERVAL"] != 0 } {incr stepnum; return} ;# skip evaluation if not on an evaluation timestep
  set incubation_time [expr $stepnum - $incubation_start]
  if { ([expr "[getstep] % $RESTART_FREQ"] == 0) && ($aborting == "False") } { ;# we need to keep track of the step number we are on
    if {$STAGE == "reverse"} {
      append_to_file $REV_FILENAME "job${ID}_${VEL_ID} $incubation_time"
    } elseif {$STAGE == "forward"} { ;# then its the forward stage
      if { $first_milestone_hit == "False" } { ;# we want to keep track of whether the starting milestone was ever crossed
        set start_crossed_string "start_uncrossed"
      } else {
        set start_crossed_string "start_crossed"
      }
      append_to_file $FWD_FILENAME "job${ID}_${VEL_ID} $incubation_time $start_crossed_string"
    } else {
      outprint "Unknown stage: $STAGE"
    }
  }
  # get coordinates of important atom groups
  loadcoords coords
  loadmasses masses
  set ligcoords [getcoords $ligatoms]
  set reccoords [getcoords $recatoms]

  # everything to perform on the first step
  if { ([getstep] == 0) || ($stepnum == 0) } {
    #set first_milestone_hit "False"
    set aborting False
    switch_grid_point $whoami
    set ligrefcoords $ligcoords
    set recrefcoords $reccoords
    set ligweights [getweights $ligatoms]
    set recweights [getweights $recatoms]
    set ligrefCOM [getCOM $ligcoords $ligweights]
    set recrefCOM [getCOM $reccoords $recweights]
    outprint "whoami: $whoami"
    outprint "MILESTONES"
    set counter 0
    foreach i $MILESTONE_LIST {
      outprint "milestone: $i"
      array set milestone $i
      if { [info exists milestone(normal)] == 1 } { set milestone(center) [ vecadd $recrefCOM [ vecscale $milestone(normal) $milestone(distance) ] ] }
      set milestone_unset [array get milestone]; array unset milestone
      #outprint "milestone_unset: $milestone_unset"
      set MILESTONE_LIST [lreplace $MILESTONE_LIST $counter $counter $milestone_unset]
      incr counter
    }
  }
  # EVERY TIMESTEP (AT INTERVAL) EXECUTES CODE BELOW

  set counter 0
  set crossed "False"
  set ligCOM [getCOM_fast $ligweights $ligcoords ]
  set recCOM [getCOM_fast $recweights $reccoords ]

  # handle each milestone
  foreach milestone_unset $MILESTONE_LIST {
    array set milestone $milestone_unset ;# first we have to unpack the milestone array
    # check milestone center_type (how the milestone moves with the receptor)
    if {(($milestone(shape) == "sphere" ) || ($milestone(shape) == "plane" ))} {
      if { $milestone(center_type) == "atom" } { ;# then the milestone is centered on an atom group
        set milestone_center_id $milestone(center_com_id)
        set milestone_center $coords($milestone_center_id)
        if {[ info exists milestone(normal) ]} {set milestone_normal $milestone(normal)} ;# only valid for planar milestones
      } elseif { $milestone(center_type) == "coord" } { ;# then the coordinate is relative to the rec COM and PA's
        # find milestone center/normal
        #outprint "coord type milestones not yet implemented"
        if { $milestone(absolute) == True } { ;# then the milestone doesn't move with the receptor
          set milestone_center $milestone(center)
          if {[ info exists milestone(normal) ]} {set milestone_normal $milestone(normal)} ;# only valid for planar milestones
        } else { ;# then the milestones DO move with the receptor
          # calculate milestone coordinates relative to the mobile receptor
          # move the milestone center relative to the rotations
          set milestone_center [rotate_by_quat $recquat $milestone(center) $recCOM $recrefCOM]
          if {[ info exists milestone(normal) ]} { ;# with planar milestones
            set milestone_normal_plus_center [rotate_by_quat $recquat [vecadd $milestone(center) $milestone(normal) ] $recCOM $recrefCOM]
            set milestone_normal [vecsub $milestone_normal_plus_center $milestone_center] ;# ... and then subtract away the new center to get the new normal
          }
        }
      }
    }

    # check the milestone shape
    set result $milestone(side) ;# by default, there should be no transition
    if { $milestone(shape) == "sphere" } {
      set result [check_spherical_milestone $ligCOM $milestone_center $milestone(radius)]
    } elseif { $milestone(shape) == "plane" } {
      set result [check_planar_milestone $ligCOM $milestone_center $milestone_normal]
    }
    # check for a crossing event
    if {[info exists milestone(siteid)]} { ;# if in this system, the milestones contain siteids
      set milestone_siteid $milestone(siteid)
    } else { # then it hasn't been implemented yet
      set milestone_siteid 0
    }
    if { ($result != $milestone(side)) && ($milestone_siteid == $SITEID) } { ;# then a crossing event has occurred
      outprint "crossing event detected for milestone: $milestone(id), result: $result, stepnum: $stepnum incubation_time: $incubation_time"
      outprint "ligand COM: $ligCOM, receptor COM: $recCOM, site COM: $milestone_center"
      #outprint "Milestones: $MILESTONE_LIST"
      set milestone(side) $result ;# set the new side that the trajectory is on
      if {$stepnum > $first_timestep_tolerance && [getstep] > $first_timestep_tolerance && ($whoami != $milestone(id) || $CARE_ABOUT_SELF == True)} { # then this is a nonzero step and the milestone sides have already been set
        set crossed $milestone(id) ;# say that we have crossed a milestone
        if { ($STAGE == "forward") && ($first_milestone_hit == "False") && ($stepnum > 0) } {
          if { ($whoami != $milestone(id)) } {
            if { $PERFECT_STARTING_DIST == "False" } { ;# unless
            # then we have a problem, because the forward stage is supposed to hit the $whoami milestone first before anything else
              outprint "FORWARD ERROR: Forward stage trajectory crossed different milestone that the starting milestone. Crossed: $milestone(id), wanted: $whoami. This trajectory must be discarded. Aborting..."
              append_to_file $FWD_FILENAME "job${ID}_${VEL_ID} failed crossed: $crossed incubation_time: $incubation_time"
              set aborting "True"
            }
          }
        }
        #if { $milestone(end) == True } { set aborting "True" } ;# if this is an ending milestone, then abort simulation
      } elseif { $whoami == $milestone(id) && $STAGE == "forward" && $first_milestone_hit == "False" && $stepnum > 0 } {
        outprint "Forward stage crossed first positional milestone. Resetting incubation time which was at: $incubation_time"
        set first_milestone_hit True
        set incubation_start $stepnum
      }
    }
    # replace the member of the list with the milestone array
    set milestone_unset [array get milestone]; array unset milestone
    set MILESTONE_LIST [lreplace $MILESTONE_LIST $counter $counter $milestone_unset]
    incr counter
  }

  if {$crossed != "False"} { ;# print useful info and, if necessary, abort simulation
    if { $aborting != "True" } {
      if { $STAGE == "reverse"} {
        set aborting "True"
        if {$whoami != $crossed } {
          outprint "Reverse succeeded."
          append_to_file $REV_FILENAME "job${ID}_${VEL_ID} succeeded incubation_time: $incubation_time"
          set reverse_failed "False"  ;# then the reversal has crossed another milestone, and we should record it
        } else {
          outprint "Reversal failed."
          append_to_file $REV_FILENAME "job${ID}_${VEL_ID} failed incubation_time: $incubation_time"
          set reverse_failed "True"
        }
      } elseif { $STAGE == "forward" } { ;# then its a forward-stage
        append_to_file $FWD_FILENAME "job${ID}_${VEL_ID} succeeded crossed: $crossed incubation_time: $incubation_time"
        outprint "GRID TRANSITION: source: $whoami, destination: $crossed, stepnum: $stepnum, ligand COM: $ligCOM, receptor COM: $recCOM, receptor start COM: $recrefCOM, incubation time: $incubation_time timesteps, incubation_start: $incubation_start, ID: $ID, VEL_ID: $VEL_ID, stage: $STAGE"

        switch_grid_point "$crossed" ;# only if its a forward stage, switch the grid point
      }
      set incubation_start $stepnum
    }
    if { ($ABORT_ON_CROSSING == True) || ($milestone(end) == True) || (($stepnum > $MAX_NUM_STEPS) && ($MAX_NUM_STEPS != 0))} { set aborting "True" }
  }
  incr stepnum
  return
}
