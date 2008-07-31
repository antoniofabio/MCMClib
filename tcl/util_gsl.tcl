proc v2l {v} {
	set ans [list]
	set d [gsl_vector_size_get $v]
	for {set i 0} {$i < $d} {incr i} {
		lappend ans [gsl_vector_get $v $i]
	}
	return $ans
}

proc l2v {l} {
	set ans [gsl_vector_alloc [llength $l]]
	for {set i 0} {$i < [llength $l]} {incr i} {
		gsl_vector_set $ans $i [lindex $l $i]
	}
	return $ans
}
