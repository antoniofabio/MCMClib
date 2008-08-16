proc stopifnot {condition} {
	if "!($condition)" {
		error "error: `$condition' is false"
	}
}

proc stopifne {a b} {
	if {[string compare $a $b]} {
		error "error: $a != $b"
	}
}

#convert a list into a double array, and back
proc ll2a {ll} {
	set ans [new_doubleArray [llength $ll]]
	for {set i 0} {$i < [llength $ll]} {incr i} {
		doubleArray_setitem $ans $i [lindex $ll $i]
	}
	return $ans
}
proc a2ll {a} {
	set ans [list]
	for {set i 0} {$i < [llength $ll]} {incr i} {
		lappend ans [doubleArray_getitem $ans $i]
	}
	return $ans
}

source ./util_gsl.tcl
source ./util_vector.tcl
