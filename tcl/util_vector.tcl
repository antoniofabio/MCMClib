#convert a list of lists into a 'vector_list' object
proc ll2vl {ll} {
	set ans [mcmclib_vector_list_alloc]
	set current $ans
	vector_list_str_v_set $current [l2v [lindex $ll 0]]
	set ll [lrange $ll 1 end]
	foreach le $ll {
		set current [mcmclib_vector_list_append [l2v $le] $current]
	}
	return $ans
}

#convert a 'vector_list' into a tcl list of lists
proc vl2ll {vl} {
	set ans [list]
	while 1 {
		lappend ans [v2l [vector_list_str_v_get $vl]]
		set vl [vector_list_str_next_get $vl]
		if {[string equal $vl NULL]} break
	}
	return $ans
}

#convert a tcl list of lists into a gsl matrix
proc ll2m {ll} {
	return [mcmclib_vector_list_asmatrix [ll2vl $ll]]
}

#convert a gsl matrix into a list of lists
proc m2ll {m} {
	set ans [list]
	set nr [gsl_matrix_size1_get $m]
	set nc [gsl_matrix_size2_get $m]
	for {set i 0} {$i < $nr} {incr i} {
		set row [list]
		for {set j 0} {$j < $nc} {incr j} {
			lappend row [gsl_matrix_get $m $i $j]
		}
		lappend ans $row
	}
	return $ans
}
