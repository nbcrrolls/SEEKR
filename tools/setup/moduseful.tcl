proc isList {arg} {	;# type checking: returns true if the argument is a list
	if { [catch {uplevel 1 "trace var $arg w {llength \[set $arg\] ;#}"}]} {
		return 1 ;# is a list
	} else {
		return 0 ;# not a list
	}
}

proc list_map {function ilist1 ilist2} {
	#this function takes two lists, and maps a function for each of them, returning the result
	#assert lengths of lists are equal
	#puts mark0
	if {[llength $ilist1] != [llength $ilist2]} {error "Error in list_map function: lists must be of equal length"}
	#type check to see whether elements of the argument lists aren't lists themselves
	#set argislist [isList [lindex $ilist1 0]]
	#puts "arg is list: $argislist"
	set bottomlevel [expr ![isList [lindex $ilist1 0]]] ;#setting the bottomlevel value to opposite the result of islist, because if islist is true, then bottomlevel needs to be false
	#puts "Bottomlevel: $bottomlevel"
	#if lowest dimension of our list
	#puts mark1
	if {$bottomlevel == 1} {
		#puts mark2.1
		set newlist {}	;#initialize the list that will contain every value within this dimension
		for {set i 0} {$i < [llength $ilist1]} {incr i} {	;#for each value in the two lists
			set result [$function [lindex $ilist1 $i] [lindex $ilist2 $i]]	;#perform the function on the two values
			lappend newlist $result 	;# append the results to be the new list's value at that location
		}
		return $newlist	;# return the dimension
	} else {	;#otherwise if we are not at the bottom level, then we'll need to recursively call the next dimension for each item
		#puts mark2.2
		set newlist {}
		for {set i 0} {$i < [llength $ilist1]} {incr i} {	;#for each descending list
			lappend newlist [list_map $function [lindex $ilist1 $i] [lindex $ilist2 $i]] 	;# recursively call this function
		}
		return $newlist ;# move up a level
	}
}

proc map {function ilist} {
	# maps a function onto a list, returning a new list with the function performed on every item in the old list
	# in order for it to work on multi-dimensional lists, need recursion for each dimension
	set bottomlevel [expr ![isList [lindex $ilist 0]]] ;#setting the bottomlevel value to opposite the result of islist, because if islist is true, then bottomlevel needs to be false
	if {$bottomlevel == 1} {
		set newlist {}
		for {set i 0} {$i < [llength $ilist]} {incr i} {	;#for each value in the list
			set result [$function [lindex $ilist $i]]	;#perform the function on the value
			lappend newlist $result 	;# append the results to be the new list's value at that location
		}
		return $newlist
	} else { ;# we need to go down another dimension
		set newlist {}
		for {set i 0} {$i < [llength $ilist]} {incr i} {	;#for each descending list
			lappend newlist [map $function [lindex $ilist $i]] 	;# recursively call this function
		}
		return $newlist ;# move up a level
	}
}

proc scientific_notation {sci_num} {
	#function converts a scientific notation of the form: NOTE: OBSELETE!!, use format command
	#  -0.000E+0.0
	#into a float number and returns it
	set numlist [split $sci_num "e|E"]
	set mantissa [lindex $numlist 0]
	set exponent [lindex $numlist 1]
	set final_result [expr "$mantissa*(pow(10,$exponent))"]
	return $final_result
}

proc heapsort_add {newnumber heap} {
	# adds an item to a heap and sorts it so that the lowest value is always in index 0
	set newheap [lappend heap $newnumber] ;# the new number goes to the end of the list
	set newheaplen [llength $newheap]
	set curindex [expr "$newheaplen - 1"]
	set nextindex [expr "($curindex - 1) / 2"]
	set nextval [lindex $newheap $nextindex]
	while {$nextval < $newnumber} {
		# we are going to swap them
		set newheap [lreplace $newheap $curindex $curindex $nextval] ;#put the upper number into the current spot
		set newheap [lreplace $newheap $nextindex $nextindex $newnumber] ;#put our new number into the higher spot
		set curindex $nextindex
		set nextindex [expr "($curindex - 1) / 2"]
		set nextval [lindex $newheap $nextindex]
	}
	return $newheap
}

proc heapsort_remove {heap} {
	# returns a heap where the lowest value at index 0 has been removed
	set endelement [lindex $heap end]
	set newheap [lreplace $heap end end] ;# delete the last item in the list
	set newheap [lreplace 0 0 $endelement] ;# replace the first item on the list with endelement
	
}

proc mergesort_by_index {seq {index 0}} {
	# recursively sorts a list of coordinates by a specific index, by devault the zeroth. OBSOLETE!! Use lsort -index instead
	set seqlen [llength $seq]
	#puts $seq
	if {$seqlen <= 1} {
		return $seq ;# its an element all by itself, it can be returned as-is
	} else { # split the list in two and recurse each element
		set half1 [lrange $seq 0 [expr "(($seqlen - 1) / 2)"]]
		set half2 [lrange $seq [expr "(($seqlen - 1) / 2) + 1"] end]
		#puts $half1
		#puts $half2
		set sortedhalf1 [mergesort_by_index $half1 $index]
		set sortedhalf2 [mergesort_by_index $half2 $index]
		#puts $sortedhalf1
		#puts $sortedhalf2
		set half1len [llength $sortedhalf1]
		set half2len [llength $sortedhalf2]
		set totallen [expr "$half1len + $half2len"]
		set totallist {}
		
		set ptr1 0
		set ptr2 0
		for {set i 0} {$i < $totallen} {incr i} {
			if {$ptr2 >= $half2len} {
				# add the remainder of sortedhalf1 to totallist
				set remainder [lrange $sortedhalf1 $ptr1 end]
				set totallist [concat $totallist $remainder]
				break
			} elseif {$ptr1 >= $half1len} {
				# add the remainder of sortedhalf2 to totallist
				set remainder [lrange $sortedhalf2 $ptr2 end]
				set totallist [concat $totallist $remainder]
				break
			} else {
			
				if {[lindex [lindex $sortedhalf1 $ptr1] $index] <= [lindex [lindex $sortedhalf2 $ptr2] $index]} { ;# then take the element from sortedhalf1
					lappend totallist [lindex $sortedhalf1 $ptr1]
					puts [lindex $sortedhalf1 $ptr1]
					incr ptr1 ;# set ptr1 to point to the next item in the list
				} else {
					lappend totallist [lindex $sortedhalf2 $ptr2]
					puts [lindex $sortedhalf2 $ptr2]
					incr ptr2 ;# set ptr2 to point to the next item in the list
				}
			}
		}
		#puts $totallist
		return $totallist ;# return the sorted list for the next round up
	}			
}








proc ClosestPair3d {pointset} {
	# given a set of 3d points, finds the closest pair
	#first construct Px and Py
	set Px [lsort -index 0 -real $pointset]
	set Py [lsort -index 1 -real $pointset]
	return [ClosestPairRec3d $Px $Py]
}

proc ClosestPairRec3d {Px, Py} {
	# recursive function for Closest pair 3d algorithm
	
}

proc slice_matrix {matrix x1 y1 x2 y2} {
  for {set i $y1} {$i < $y2} {incr i} {
    lappend matrixslice [lrange [lindex $matrix $i] $x1 [expr "$x2 - 1"]] ;# get region A and put it into the oldmatrix
  }
  return $matrixslice
}

proc join_matrix {matrix1 matrix2} {
  if {[llength $matrix1] != [llength $matrix2]} {error "join_matrix: Matrices are not same height"} ;# assertion
  for {set i 0} {$i < [llength $matrix1]} {incr i} {
    lappend newmatrix [concat [lindex $matrix1 $i] [lindex $matrix2 $i]]
  }
  return $newmatrix
}

proc writepqr {sel filename} { ;# writes a pqr file for input into Delphi
	#puts mark1
	# gets the lists about the selection
	set record_name "ATOM  "
	set atom_number_list [$sel get "serial"] ;# because serial starts with number 1
	set atom_name_list [$sel get "name"]
	set residue_name_list [$sel get "resname"]
	set chain_id_list [$sel get "chain"]
	set residue_number_list [$sel get "resid"]
	#puts mark5
	set space "    "
	set xvalue_list [$sel get "x"]
	set yvalue_list [$sel get "y"]
	set zvalue_list [$sel get "z"]
	set charge_list [$sel get "charge"]
	set radius_list [$sel get "occupancy"]
	set file_data {}
	# gets all parameters from the selected atom
	#puts mark9
	set atomcount [llength $atom_number_list]
	for {set i 0} {$i < $atomcount} {incr i} {
		set atom_number [lindex $atom_number_list $i]
		set atom_name [lindex $atom_name_list $i]
		set residue_name [lindex $residue_name_list $i]
		set chain_id [lindex $chain_id_list $i]
		set residue_number [lindex $residue_number_list $i]
		set xvalue [lindex $xvalue_list $i]
		set yvalue [lindex $yvalue_list $i]
		set zvalue [lindex $zvalue_list $i]
		set charge [lindex $charge_list $i]
		set radius [lindex $radius_list $i]
		# writes the atoms in pqr format
		set newline [format "%6s%5d%5s%4s%2s%4d%4s%8.3f%8.3f%8.3f%8.4f%7.4f\n" $record_name $atom_number $atom_name $residue_name $chain_id $residue_number $space $xvalue $yvalue $zvalue $charge $radius]
		lappend filedata $newline
	}
	set filedata [join $filedata ""]
	# open the file and everything
	set pqrfile [open $filename w]
	puts $pqrfile "REMARK   1 PQR file generated by DelEnsembleElec"
	puts $pqrfile $filedata
	close $pqrfile
	
}

proc set_append {ilist appendee} {
  #returns ilist with appendee as a member only if appendee didn't already exist there, like a mathematical set
  if {[lsearch $ilist $appendee] == -1} {
    #it doesn't exist in the set
    lappend ilist $appendee
  }
  #could try converting into a string to remove those weird marks
  return $ilist
}

proc mean {ilist} {
  # returns a value which is the mean of every element in the list
  # for map functions
  set sum 0
  set listlen [llength $ilist]
  foreach member $ilist {
    set sum [expr "$sum + $member"]
  }
  set avg [expr "$sum / $listlen"]
  return $avg
}

proc slice_matrix {matrix x1 y1 x2 y2} {
  #slices and returns a rectangle from a 2d array
  set matrixslice {}
  if {$x2 == "end"} {set xlast "end"} else {set xlast [expr "$x2 - 1"]} 
  if {$y2 == "end"} {set y2 [llength $matrix]}
  for {set i $y1} {$i < $y2} {incr i} {
    lappend matrixslice [lrange [lindex $matrix $i] $x1 $xlast] ;# get region A and put it into the oldmatrix
  }
  return $matrixslice
}

proc join_matrix {matrix1 matrix2} {
  set newmatrix {}
  if {[llength $matrix1] != [llength $matrix2]} {error "join_matrix: Matrices are not same height"} ;# assertion
  for {set i 0} {$i < [llength $matrix1]} {incr i} {
    lappend newmatrix [concat [lindex $matrix1 $i] [lindex $matrix2 $i]]
  }
  return $newmatrix
}

proc gaussian3d {x y z x0 y0 z0 sigma A} {
	# NOTE: sigma is considered the same in all 3 dimensions
	set sigsquared2 [expr "2.0*($sigma**2)"]
	set xterm [expr "(($x - $x0)**2)/($sigsquared2)"]
	set yterm [expr "(($y - $y0)**2)/($sigsquared2)"]
	set zterm [expr "(($z - $z0)**2)/($sigsquared2)"]
	set answer [expr "$A * ($e**-($xterm + $yterm + $zterm))"]
	return $answer	
}

proc interp_sphere {four_coords} {
  # given 4 coordinates, will find the equation for a sphere that interpolates all of them
  if {[llength $four_coords] != '4'} {error "interp_sphere argument four_coords but be a list of length 4"}
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
 
proc LU_factorize {A N} {
  if {$N == 1} { ;# then matrix A is 1x1, we can return it
    return $A
  }
  puts mark0
  set top_pivot [lindex [lindex $A 0] 0] ;# top leftmost member of matrix
  puts mark1
  set inv_top_pivot [expr "1 / $top_pivot"] ;# a little speedup by avoiding dividing on iterations
  puts mark2
  for {set row 1} {$row < $N} {incr row} { ;# for each row except the top
    puts "A: $A"
    set pivot [lindex [lindex $A $row] 0] ;# leftmost member of matrix
    set multiplier [expr "$pivot * $inv_top_pivot"]
    puts mark3
    set U_little_row {}
    for {set col 1} {$col < $N} {incr col} {
      # multiply first row times multiplier, then subtract from this row
      lappend U_little_row [expr "[ lindex [ lindex $A $row] $col] - ([ lindex [ lindex $A 0] $col] * $multiplier)"]
      puts mark4
      puts "U: $U_little_row"
    }
    lappend U_little $U_little_row ;# add the row to the U_little matrix
    lappend multipliers $multiplier
  }
  puts mark5
  set U_little_solved [LU_factorize $U_little [expr "$N - 1"]] ;# recursively solve the sub-matrix
  lappend U [lindex $A 0] ;# append the top row of A to U
  for {set row 0} {$row < [expr "$N - 1"]} {incr row} {
    lappend U "[lindex $multipliers $row] [lindex $U_little_solved $row]"
  }
  return $U ;# return factorized matrix    
}

proc test_time {command} {
  set starttime [clock microseconds]
  eval $command
  return [ expr "[clock microseconds] - $starttime"]
}

proc householder { matrix } {
  # performs the householder transformation on a symmetric matrix
  set n [llength $matrix]
  for {set k 0} {$k < [expr "$n - 2"]} {incr k} {
    # calculate q
    set q 0.0; set v ""; set u ""; set z ""
    for {set j [expr "$k + 1"]} {$j < $n} {incr j} {
      set a [lindex [lindex $matrix $j] $k]
      set q [expr "$q +  $a*$a "]
    }
    # calculate alpha
    set a [lindex [lindex $matrix [expr "$k + 1"]] $k]
    if { $a==0 } {
      set alpha [expr "-sqrt($q)"]
    } else {
      set alpha [expr "-sqrt($q)*$a / abs($a)"]
    }
    set RSQ [expr "$alpha * $alpha - $alpha * $a"] ;# RSQ = 2*r**2
    set r [expr "sqrt($RSQ/2.0)"]
    # populate v vector
    lappend v 0
    lappend v [expr "$a - $alpha"]
    set w [expr "($a - $alpha) / (2.0 * $r)"]; puts "w1: $w"
    for {set j [expr "$k + 2"]} {$j < $n} {incr j} {
      set a [lindex [lindex $matrix $j] $k]
      lappend v $a
      set w [expr "$a / (2.0 * $r)"]
    }
    for {set j $k} {$j < $n} {incr j} {
      set j_index [expr "$j - $k"]
    }
    # populate u vector
    for {set j $k} {$j < $n} {incr j} {
      set uj 0
      for {set i [expr "$k+1"]} {$i < $n} {incr i} { ;# summation
        set vi [lindex $v [expr "$i - ($k)"]]
        set uj [expr "$uj + [lindex [lindex $matrix $j] $i] * $vi"]
      }
      set uj [expr "$uj / $RSQ"]
      lappend u $uj
    }
    # get dot product of the two vectors
    set PROD [cheap_vecdot $u $v]
    # populate z vector
    for {set j $k} {$j < $n} {incr j} {
      set uj [lindex $u [expr "$j - $k"]]
      set vj [lindex $v [expr "$j - $k"]]
      set zj [expr "$uj - ($PROD/(2*$RSQ))*$vj"]
      lappend z $zj
    }
    # replace values in matrix
    for { set ell [expr "$k + 1"] } {$ell < [expr "$n - 1"]} {incr ell} { ;# I'm using 'ell' because I don't want to use merely the letter 'l'; looks too much like a '1'
      set ell_index [expr "$ell - ($k)"]
      for {set j [expr "$ell + 1"]} {$j < $n} {incr j} {
        set j_index [expr "$j - ($k)"]
        set value [expr "[lindex [lindex $matrix $j] $ell] - [lindex $v $ell_index]*[lindex $z $j_index] - [lindex $v $j_index]*[lindex $z $ell_index]"]
        lset matrix $ell $j $value; lset matrix $j $ell $value
      }
      set value [expr "[lindex [lindex $matrix $ell] $ell] - 2.0 * [lindex $v $ell_index]*[lindex $z $ell_index]"]
      lset matrix $ell $ell $value
    }
    # replace other value in matrix
    set n_1 [expr "$n - 1"]
    lset matrix $n_1 $n_1 [expr "[lindex [lindex $matrix $n_1] $n_1] - 2.0 * [lindex $v end]*[lindex $z end]"]
    for {set j [expr "$k+2"]} {$j < $n} {incr j} {
      lset matrix $k $j 0; lset matrix $j $k 0
    }
    set k_plus [expr "$k+1"]
    set value [expr "[lindex [lindex $matrix $k] $k_plus] - [lindex $v 1]*[lindex $z 0]"]
    lset matrix $k $k_plus $value; lset matrix $k_plus $k $value
  }
  return $matrix
}

proc largest_eigenvector_QR { matrix { TOL 1.0e-4 } { MAXITER 1000 } } {
  # Uses QR factorization to find the largest eigenvectors/eigenvalues of a symmetric tridiagonal matrix
  set n [llength $matrix]
  set a "[lindex [lindex $matrix 0] 0]"; set b "0.0" ;# we don't need the whole matrix, just the main diagonal and the first subdiagonal
  set evals []; set evecs []
  for {set i 1} {$i < $n} {incr i} {
    lappend a [lindex [lindex $matrix $i] $i]
    lappend b [lindex [lindex $matrix $i] [expr "$i - 1"]]
    puts "a$i: [lindex $a end] b$i: [lindex $b end]"
  }
  set SHIFT 0
  for {set k 1} {$k < $MAXITER} {incr k} {
    if {[expr "abs([lindex $b end])"] <= $TOL} { ;# check the bottom corner
      set lambda [expr "[lindex $a end] + $SHIFT"]
      lappend evals $lambda
      set n [expr "$n - 1"]
    }
    if {[expr "abs([lindex $b 1])"] <= $TOL} { ;# check the top corner
      set lambda [expr "[lindex $a 0] + $SHIFT"]
      lappend evals $lambda
      set n [expr "$n - 1"]
      lset a 0 [lindex $a 1]
      for {set j 1} {$j < $n} {incr j} {
        lset a $j [lindex $a [expr "$j + 1"]]
        lset b $j [lindex $b [expr "$j + 1"]]
      }
    }
    if {$n == 0} {
      return "$evals $evecs"
    }
    for {set j 2} {$j < [expr "$n - 1"]} {incr j} {
      if {[expr "abs([lindex $b $j])"] <= $TOL} {
        
      }
    }
    #NOTE: INCOMPLETE!!!!!!!!!!!
  }
}

proc largest_eigenvector_powermethod { matrix {guess {1 0 0 0}} {TOL 1.0e-4} {MAXITER 10000}} {
  # calculates the principle eigenvector of matrix using the Symmetric Power Method
  set x [cheap_vecnorm $guess] ;# normalize the guess vector
  #set mu0 0 ;# used for Aitken's method to speed convergence; All reference of mu commented because those are e.values and unnecessary for the calculation
  #set mu1 0
  set counter 0
  for {set k 0} {$k < $MAXITER} {incr k} {
    set y [mat_vec_mult $matrix $x]
    #set mu [cheap_vecdot $x $y] ;# approximation of the eigenvalue
    #puts "mu: $mu, mu0: $mu0, mu1: $mu1"
    #if {($mu != 0) || ($mu0 != 0) || ($mu1 != 0)} {set mu_hat [expr "$mu0 - (($mu1 - $mu0)*($mu1 - $mu0)/($mu - 2*$mu1 + $mu0))"]}
    if {[veclength $y] == 0} {error "cannot find reasonable eigenvector in function 'largest_eigenvector', please try different guess and resubmit"}
    set normed_y [cheap_vecnorm $y]
    set err [veclength [vecsub $x $normed_y]]
    set x $normed_y
    #puts "x: $x, err: $err"
    if {$err < $TOL} {puts "iter: $counter"; return $x} ;#return the closest eigenvalue
    #set mu0 $mu1; set mu1 $mu
    incr counter
  }
  error "function: 'largest_eigenvector' could not converge to an eigenvalue within the iterations specified for matrix: $matrix"
}

proc eye_vec {coords {molid 0}} {
  # coordinates of the atom
  #set coords [lindex [$sel get {x y z}] 0]
  set ourmol $molid ;#[$sel molid]
  # position in world space
  set mat [lindex [molinfo $ourmol get view_matrix] 0]
  set world [vectrans $mat $coords]

  # since this is orthographic, I just get the projection
  lassign $world x y
  # get a coordinate behind the eye
  set world2 "$x $y 5"

  # convert back to molecule space
  # (need an inverse, which is only available with the measure command)
  set inv [measure inverse $mat]
  set coords2 [vectrans $inv $world2]
  draw cylinder $coords $coords2 radius 0.3
  return [vecsub $coords2 $coords]
}

proc draw_circle {center radius {resolution 100}} {
  set PI 3.141592
  set eyevec [eye_vec $center]
  set ortho [vecscale [vecnorm [veccross $eyevec "1 0 0"]] $radius]
  set dtheta [expr "2 * $PI / $resolution"]
  set rotmat [trans center $center axis $eyevec $dtheta rad]
  set oldortho $ortho
  for {set i 0} {$i < $resolution} {incr i} {
    set ortho [vectrans $rotmat $ortho]
    draw line [vecadd $center $ortho] [vecadd $center $oldortho]
    
    set oldortho $ortho
  }
  return 
}
