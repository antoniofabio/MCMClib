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
