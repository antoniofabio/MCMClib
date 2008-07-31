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
