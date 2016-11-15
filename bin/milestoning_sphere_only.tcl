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

if {!([info exists ID])} {set ID 0}
if {!([info exists SITEID])} {set SITEID 0}
#set LA_src "$TEMPLATE_LA_src" ;#NOTE: this will need to be defined by the user
#if {$LA_src != ""} { source $LA_src } ;# import the tcl linear algebra package

set incubation_start 0


proc add_id_to_FHPD {id} {
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

proc closest_grid_point {ligcenter neighbors} {
  ###########################################################################################
  # finds point in $neighbors closest to $ligcenter.
  ###########################################################################################
  array set array_neighbors $neighbors ;# pass argument as a list, convert to array once we're here
  set closestdist 999999
  set closestneighbor None
  foreach neighbor [array names array_neighbors] { ;# loop through each neighbor to see which one we are closest to
    #puts "neighbor: $neighbor, $array_neighbors($neighbor)"
    set neighbor_lig_dist [cheap_vecdist $neighbor $ligcenter]
    if {$neighbor_lig_dist < $closestdist} { ;# set this to be the closest neighbor
      set closestdist $neighbor_lig_dist
      set closestneighbor $array_neighbors($neighbor)
    }
  }
  return $closestneighbor
}

proc orient {center vector1 vector2} { ;# given a center and 2 vectors, will return the matrix to align vector 1 to vector 2
  set vec1 [cheap_vecnorm $vector1]
  set vec2 [cheap_vecnorm $vector2]
  set rotvec [veccross $vec1 $vec2]
  set sine   [veclength $rotvec]
  set cosine [vecdot $vec1 $vec2]
  set angle [expr atan2($sine,$cosine)]  
  # return the rotation matrix
  set mat [cheap_trans $center $rotvec $angle]
  return $mat
}

proc easy_principalaxes {COM PA1_group PA3_group} {
  # Given two groups whose COM roughly lies on the principal axes, function will
  # find pseudo-PA with 'easier' or quicker calculations
  #output "PA1_group: $PA1_group"
  #output "COM: $COM"
  set PA1 [cheap_vecnorm [vecsub $PA1_group $COM]] ;# this will be the first PA
  set PA3_skew [vecsub $PA3_group $COM] ;# this is just for the cross product
  set PA2 [cheap_vecnorm [veccross $PA1 $PA3_skew]] ;# 2nd PA
  set PA3 [cheap_vecnorm [veccross $PA1 $PA2]] ;# 3rd PA
  return "{$PA1} {$PA2} {$PA3}"
}

proc planar_regression {allcoords {dep_axiz z}} {
  # given an atomselection, will calculate a planar regression
  set sumx2 0.0; set sumy2 0.0; set sumxy 0.0; set sumxz 0.0
  set sumyz 0.0; set sumx 0.0; set sumy 0.0; set sumz 0.0
  set n [llength $allcoords]
  foreach coord $allcoords {
    set x [lindex $coord 0]
    set y [lindex $coord 1]
    set z [lindex $coord 2]
    set sumx2 [expr "$sumx2 + ($x * $x)"]
    set sumy2 [expr "$sumy2 + ($y * $y)"]
    set sumxy [expr "$sumxy + ($x * $y)"]
    set sumxz [expr "$sumxz + ($x * $z)"]
    set sumyz [expr "$sumyz + ($y * $z)"]
    set sumx [expr "$sumx + $x"]
    set sumy [expr "$sumy + $y"]
    set sumz [expr "$sumz + $z"]
  }
  # matrix A: 3x3 matrix
  set A "2 3 3 $sumx2 $sumxy $sumx $sumxy $sumy2 $sumy $sumx $sumy $n"
  set b "2 3 1 $sumxz $sumyz $sumz"
  set soln [La::msolve $A $b]
  set eq [lrange $soln 3 end]
  set PA1 [cheap_vecnorm "1 0 [lindex $eq 0]"]
  set PA2 [cheap_vecnorm "0 1 [lindex $eq 1]"]
  set PA3 [cheap_vecnorm [veccross $PA1 $PA2]]
  return "{$PA1} {$PA2} {$PA3}"
}

proc output { message } { ;# automatically appends caption to stdout from script
  global outputprefix;
  puts "$outputprefix$message"
}

proc range {start cutoff finish {step 1}} {
	# If "start" and "finish" aren't integers, do nothing:
	if {[string is integer -strict $start] == 0 || [string is\
		integer -strict $finish] == 0} {
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
		        error "range: Use \"to\" for an inclusive\
		        range, or \"no\" for a noninclusive range"
		}
	}

	# Is the range ascending or descending (or neither)?
	set up [expr {$start <= $finish}]

	# If range is descending and step is positive but
	# doesn't have a "+" sign, change step to negative:
	if {$up == 0 && $step > 0 && [string first "+" $start] != 0} {
		set step [expr $step * -1]
	}

	# Initialize list variable and identify
	# class of integer range:
	set ranger [list]
	switch "$up $inclu" {
		"1 1" {set op "<=" ; # Ascending, inclusive range}
		"1 0" {set op "<" ;  # Ascending, noninclusive range}
		"0 1" {set op ">=" ;  # Descending, inclusive range}
		"0 0" {set op ">" ;  # Descending, noninclusive range}
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

proc cheap_veccross {a b} { ;# cross product between two vectors (of length 3) DEPRECATED
  set a1 [lindex $a 0]; set a2 [lindex $a 1]; set a3 [lindex $a 2]
  set b1 [lindex $b 0]; set b2 [lindex $b 1]; set b3 [lindex $b 2]
  set v1 [expr "$a2*$b3 - $a3*$b2"]
  set v2 [expr "$a3*$b1 - $a1*$b3"]
  set v3 [expr "$a1*$b2 - $a2*$b1"]
  return "$v1 $v2 $v3"
}

proc cheap_vecnorm {a} { ;# normalize a vector to be of length 1.0
  set scalefactor [expr "1/[veclength $a]"]
  return [vecscale $a $scalefactor]
}

proc cheap_vecdot {a b} { ;# dot two vectors onto each other DEPRECATED
  set veclen [llength $a]
  if {$veclen != [llength $b]} {error "cheap_vecdot: Dot product attempted between two vectors of different lengths: $a dot $b"}
  set total 0.0
  for {set i 0} {$i < $veclen} {incr i} {
    set total [expr "$total + [lindex $a $i]*[lindex $b $i]"]
  }
  return $total
}

proc mat_mat_mult {m1 m2} { ;# performs matrix matrix multiplication m1 * m2 DEPRECATED
  set m1width [llength [lindex $m1 1]]
  set m1height [llength $m1]
  set m2width [llength [lindex $m2 1]]
  set m2height [llength $m2]
  if {$m1width != $m2height} {error "matrices not aligned for matrix-matrix multiplication"}
  set m { }
  for {set i 0} {$i < $m1height} {incr i} {
    set row { }
    for {set j 0} {$j < $m2width} {incr j} { ;# for every entry in the resulting matrix
      set entry 0.0
      for {set k 0} {$k < $m1width} {incr k} { ;# for every entry in the source matrices
        set entry [expr "$entry + [lindex [lindex $m1 $i] $k]*[lindex [lindex $m2 $k] $j]"]
      }
      lappend row $entry
    }
    lappend m $row
  }
  return $m
}

proc mat_vec_mult {mat vec} { ;#DEPRECATED
  set mwidth [llength [lindex $mat 1]]
  set mheight [llength $mat]
  set vlen [llength $vec]
  set newvec { }
  if {$mwidth != $vlen} {error "matrix/vector not properly aligned for multiplication"}
  for {set i 0} {$i < $mheight} {incr i} {
    set entry 0.0
    for {set j 0} {$j < $vlen} {incr j} {
      set entry [expr "$entry + [lindex [lindex $mat $i] $j]*[lindex $vec $j]"]
    }
    lappend newvec $entry
  }
  return $newvec
}

proc translate_mat {transvec} {
  # create a matrix to translate a point somewhere else
  return "{1.0 0 0 [lindex $transvec 0]} {0 1.0 0 [lindex $transvec 1]} {0 0 1.0 [lindex $transvec 2]} {0 0 0 1.0}"
}

proc cheap_trans {center rotvec angle} { ;# given a center vector, rotation
  #centering
  set t1 "{1.0 0 0 -[lindex $center 0]} {0 1.0 0 -[lindex $center 1]} {0 0 1.0 -[lindex $center 2]} {0 0 0 1.0}"
  #rotation
  #puts "veclength: [veclength $rotvec]"
  if {[veclength $rotvec] == 0.0} {set rotvec {1 1 1}; set angle 0.0} ;# arbitrary rotation vector and no angle: no transform
  set u [cheap_vecnorm $rotvec] ;# normalize the rotation vector
  set ux [lindex $u 0]; set uy [lindex $u 1]; set uz [lindex $u 2]
  set sintheta [expr "sin($angle)"]
  set costheta [expr "cos($angle)"]
  set one_minus_costheta [expr "1-$costheta"]
  set r0_0 [expr "$costheta + $ux*$ux*$one_minus_costheta"]
  set r0_1 [expr "$ux*$uy*$one_minus_costheta - $uz*$sintheta"]
  set r0_2 [expr "$ux*$uz*$one_minus_costheta + $uy*$sintheta"]
  set r1_0 [expr "$uy*$ux*$one_minus_costheta + $uz*$sintheta"] ;# filling out the rotation portion of the matrix
  set r1_1 [expr "$costheta + $uy*$uy*$one_minus_costheta"]
  set r1_2 [expr "$uy*$uz*$one_minus_costheta - $ux*$sintheta"]
  set r2_0 [expr "$uz*$ux*$one_minus_costheta - $uy*$sintheta"]
  set r2_1 [expr "$uz*$uy*$one_minus_costheta + $ux*$sintheta"]
  set r2_2 [expr "$costheta + $uz*$uz*$one_minus_costheta"]
  set r " {$r0_0 $r0_1 $r0_2 0.0} {$r1_0 $r1_1 $r1_2 0.0} {$r2_0 $r2_1 $r2_2 0.0} {0.0 0.0 0.0 1.0}"
  # reverse the centering
  set t2 "{1.0 0 0 [lindex $center 0]} {0 1.0 0 [lindex $center 1]} {0 0 1.0 [lindex $center 2]} {0 0 0 1.0}"
  # combine matrices
  set t1_r [transmult $t2 $r]
  set m [transmult $t1_r $t1]
  return $m
}

proc cheap_vecbackscale {src vec scale} { ;# given a vector placed at a source, will scale it backwards such that the source changes
  set tip [vecadd $src $vec] ;# the tip of the vector
  set backvec [vecsub $tip $src] ;# vec pointing backwards
  set newbackvec [vecscale $backvec $scale] ;# scale by the amount we want
  set newsrc [vecsub $tip $newbackvec]
  return $newsrc
}

###########################################################################################
# INITIALIZATIONS: continued
###########################################################################################

#if {[llength $LIGRANGE] > 2} {eval "set liglist \[range $LIGRANGE\]"} else {set liglist $LIGRANGE}
#if {[llength $RECRANGE] > 2} {eval "set reclist \[range $RECRANGE\]"} else {set reclist $RECRANGE} ;# we have to do all this jury-rigging with 'eval' because namd tcl won't take {*} expansion
set liglist [parse_selection $LIGRANGE]
set reclist [parse_selection $RECRANGE]
;# get center of mass of receptor/ligand
#puts "liglist: $liglist"
set ligCOMid [addgroup $liglist]
set recCOMid [addgroup $reclist]

set MILESTONE_LIST_COPY $MILESTONE_LIST
set MILESTONE_LIST ""
foreach milestone_unset $MILESTONE_LIST_COPY {
  array set milestone $milestone_unset
  if {$milestone(center_type) == "atom"} { ;# then we will want to add these atoms to a coordinate group
    set group [addgroup $milestone(center_indeces)] ;# define the group that represents the center of this milestone
    set milestone(center_com_id) $group
    
  }
  
  set milestone_unset [array get milestone]
  lappend MILESTONE_LIST $milestone_unset
}

set recidlist {}
foreach atom $reclist {
  lappend recidlist [addgroup $atom]
}
set PA1idlist {}
foreach atom $RECPA1_LIST {
  lappend PA1idlist [addgroup $atom]
}
set recPA1id [addgroup $RECPA1_LIST]
set PA3idlist {}
foreach atom $RECPA3_LIST {
  lappend PA3idlist [addgroup $atom]
}
set recPA3id [addgroup $RECPA3_LIST]

set startrecPA {}
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

proc calc_ref_movements {coord refCOM refCOMstart refPA refPAstart} {
  ###########################################################################################
  # calculates the location of a point relative to the movements of a frame of reference
  #  coord: the point relative to which the frame of reference has moved
  #  refCOM: the center of mass of the frame of reference 
  #  refCOMstart: the starting location of the frame's center of mass
  #  refPA: the principal axes of the frame
  #  refPAstart
  ###########################################################################################
  # frame of reference translation = refCOM - refCOMstart
  # frame of reference rotation = refPAstart --> refPA
  
  # compute the angle and axis of rotation
  if { $refPAstart == ""} { # then there is no reference rotation
    set refPAstart { {1 0 0} {0 1 0} {0 0 1} }
    set refPA { {1 0 0} {0 1 0} {0 0 1} }
  }
  set matrix1 [orient {0 0 0} [lindex $refPAstart 2] [lindex $refPA 2]] ;# orient z to the 3rd principal axis
  set newaxes {} ;# need to define a new set of 3 axes based on the rotation applied by matrix 1
  foreach start_axis $refPAstart {lappend newaxes [lrange [vectrans $matrix1 "$start_axis 1"] 0 2]} ;# apply matrix 1 to the starting axes
  set matrix2 [orient {0 0 0} [lindex $newaxes 1] [lindex $refPA 1]] ;# orient the new axes to the 2nd principal axis
  set translation_matrix [translate_mat [vecsub $refCOM $refCOMstart]]
  set final_matrix [transmult $translation_matrix [transmult $matrix2 $matrix1]] ;# this is the matrix needed to align the starting PA's to the current PA's  
  # reset coord to be aligned to the reference
  set newcoord [lrange [vectrans $final_matrix "$coord 1"] 0 2] ;# annoyingly, we have to add a 1 to the end of this variable for multiplcation, then remove it again afterwards
  return $newcoord
}

proc make_reversals { flag } {
  # This function checks if velocity reversals have been properly implemented in this NAMD.
  # If so, it flags NAMD to reverse velocities in the current step
  output "REVERSING ALL VELOCITIES!!!"
  if {[ catch {setboundaryflag $flag} ] != 0} {output "ERROR: setboundaryflag command not recognized in this compilation of NAMD."; abort}
}

proc calcforces { } { 
  global ligCOMid; global recCOMid; global recPA1id; global recPA3id; global recidlist; global PA1idlist; global PA3idlist; global neighbors; global raw_neighbors; global startrecPA; global startrecCOM; global gridcenter; global whoami; global incubation_start; global recCOM; global SCRIPT_INTERVAL; global PHASE; global CARE_ABOUT_SELF; global LIGROT; global ABORT_ON_CROSSING; global LIGRANGE; global RECRANGE; global RECPA1_LIST; global RECPA3_LIST; global FHPD_FILE; global MILESTONE_LIST; global first_timestep_tolerance; global ID; global first_milestone_hit; global SITEID;
  #output "TEST MARK 0"
  set stepnum [getstep] ;# step number
  setboundaryflag 0 ;# make sure the reversals are turned on only for the steps we want them turned on for
  #make_reversals 0
  if { [expr "$stepnum % $SCRIPT_INTERVAL"] != 0 } {return} ;# skip evaluation if not on an evaluation timestep
  #output "TESTING OUTPUT"
  # get coordinates of important atom groups
  loadcoords coords
  set ligCOM $coords($ligCOMid)
  set recCOM $coords($recCOMid)
  set reccoordlist {}
  foreach recid $recidlist {
    lappend reccoordlist $coords($recid)
  }
  
  #PA 1
  set recPA1coordlist {}
  foreach PA1id $PA1idlist {
    lappend recPA1coordlist $coords($PA1id)
  }
  set recPA1 $coords($recPA1id)
  
  #PA 3
  set recPA3 $coords($recPA3id) 
  set recPA3coordlist {}
  foreach PA3id $PA3idlist {
    lappend recPA3coordlist $coords($PA3id)
  }
  set recPA3 $coords($recPA3id) 
  
  # everything to perform on the first step
  if { $stepnum == 0 } { ;
    if {$PHASE == "reverse"} { make_reversals 1 } ;# for the reversal phase of the milestoning: reverse all velocities to start with. NOTE: requires special build of NAMD
    switch_grid_point $whoami
    #puts "first gridcenter: 
    if { $RECPA1_LIST != "" } {
      set startrecPA [easy_principalaxes $recCOM $recPA1 $recPA3] ;# needed for alignment of rec in each frame to the first frame
    } else {
      set startrecPA ""
    }
    set startrecCOM $recCOM
    
    output "whoami: $whoami"
    output "MILESTONES"
    set counter 0
    foreach i $MILESTONE_LIST {
      output "milestone: $i"
      array set milestone $i
      #output "startrecCOM: $startrecCOM"
      if {([ info exists milestone(normal) ]) && ([ info exists milestone(distance) ]) && ([ info exists milestone(center) ])} {set milestone(center) [ vecadd $startrecCOM [ vecscale $milestone(normal) $milestone(distance) ] ]}
      set milestone_unset [array get milestone]
      #output "milestone_unset: $milestone_unset"
      set MILESTONE_LIST [lreplace $MILESTONE_LIST $counter $counter $milestone_unset] 
      incr counter
    }
    #for {set i 0} {$i < [expr "[llength $raw_neighbors]/2"]} {incr i} {lappend neighbors [vecadd [lindex $raw_neighbors [expr "$i*2"]] $gridcenter]; lappend neighbors [lindex $raw_neighbors [expr "$i*2 + 1"]]}
  } ;# vector from rec COM with rec PAs as unit vector describing the center of the grid
  
  # handle each milestone
  set counter 0
  set crossed False; set aborting False
  foreach milestone_unset $MILESTONE_LIST {
    array set milestone $milestone_unset ;# first we have to unpack the milestone array
    #if { (($CARE_ABOUT_SELF == False) || ($stepnum <= $first_timestep_tolerance)) && ($whoami == $milestone(id))} { incr counter; continue } ;# if milestones are ignoring self-crossings, then skip this milestone. Also skip if the milestone is crossing itself, and the trajectory needed some time to move to a side. Increment the counter anyways
    
    # check milestone center_type (how the milestone moves with the receptor)
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
        if { $RECPA1_LIST != "" } {
          set recPA [easy_principalaxes $recCOM $recPA1 $recPA3] ;# needed for alignment of rec in each frame to the first frame
        } else {
          set recPA ""
        } 
        set milestone_center [ calc_ref_movements $milestone(center) $recCOM $startrecCOM $recPA $startrecPA ]
        #output "recCOM: $recCOM, recPA: $recPA, milestone_center: $milestone_center, milestone(center): $milestone(center)"
        if {[ info exists milestone(normal) ]} { ;# with planar milestones
          set milestone_normal_plus_center [ calc_ref_movements [vecadd $milestone(center) $milestone(normal)] $recCOM $startrecCOM $recPA $startrecPA ] ;# in order to get the normal, need to add normal to center, calculate shift in frame of ref.
          set milestone_normal [vecsub $milestone_normal_plus_center $milestone_center] ;# ... and then subtract away the new center to get the new normal
        }
      } 
    }
    # check the milestone shape
    #output "milestone: $milestone_unset"
    if { $milestone(shape) == "sphere" } {
      #puts "checking spherical milestone: ligCOM: $ligCOM, milestone_center: $milestone_center, milestone(radius): $milestone(radius)"
      set result [check_spherical_milestone $ligCOM $milestone_center $milestone(radius)]
    } elseif { $milestone(shape) == "plane" } {
      set result [check_planar_milestone $ligCOM $milestone_center $milestone_normal]
    } ;# elseif other milestone types ???
    # check for a crossing event
    if {[info exists milestone(siteid)]} { ;# if in this system, the milestones contain siteids
      set milestone_siteid $milestone(siteid)
    } else { # then it hasn't been implemented yet
      set milestone_siteid 0
    }
    if { ($result != $milestone(side)) && ($milestone_siteid == $SITEID) } { ;# then a crossing event has occurred
      set milestone(side) $result ;# set the new side that the trajectory is on
      if {$stepnum > $first_timestep_tolerance && ($whoami != $milestone(id) || $CARE_ABOUT_SELF == True)} { # then this is a nonzero step and the milestone sides have already been set
        set crossed $milestone(id) ;# say that we have crossed a milestone
        if { $milestone(end) == True } { set aborting True } ;# if this is an ending milestone, then abort simulation
      }
      if { ($PHASE == "forward") && ($first_milestone_hit == "False") && ($stepnum > 0) } {
        if { ($whoami != $milestone(id)) } {
          # then we have a problem, because the forward phase is supposed to hit the $whoami milestone first before anything else
          output "FORWARD ERROR: Forward phase trajectory crossed different milestone that the starting milestone. Wanted: $milestone(id), crossed: $whoami. This trajectory must be discarded. Aborting..."
          abort
        } else {
          # then the forward phase has crossed the starting plane, as expected.
          output "Forward phase crossed first milestone. Resetting incubation time which was at: [expr $stepnum - $incubation_start]"
          set first_milestone_hit True
          set incubation_start $stepnum
        }
      }
    }
    # replace the member of the list with the milestone array
    set milestone_unset [array get milestone]
    #output "milestone_unset: $milestone_unset"
    set MILESTONE_LIST [lreplace $MILESTONE_LIST $counter $counter $milestone_unset] 
    incr counter
  }
  
  #output "Milestones: $MILESTONE_LIST"
  #set crossed [barrier_check $neighbors $ligCOM $recCOM $recPA ] ;# call to check if planar/rotational barriers crossed
  #puts "made it this far"
  if {$crossed != "False"} { ;# print useful info and, if necessary, abort simulation
    output "Milestones: $MILESTONE_LIST"
    output "GRID TRANSITION: source: $whoami, destination: $crossed, current step: $stepnum, ligand COM: $ligCOM, receptor COM: $recCOM, receptor start COM: $startrecCOM, incubation time: [expr {$stepnum - $incubation_start}] timesteps"
    if { $PHASE == "reverse" && $whoami != $crossed } { output "adding to FHPD file"; add_id_to_FHPD $ID } ;# then the reversal has crossed another milestone, and we should record it
    set incubation_start $stepnum
    switch_grid_point "$crossed"
    if { ($ABORT_ON_CROSSING == True) || ($aborting == True) || (($stepnum > $MAX_NUM_STEPS) && ($MAX_NUM_STEPS != 0))} { abort }
  }
}

###########################################################################################
# optimization section, variables that have been set to optimize simulation
###########################################################################################


