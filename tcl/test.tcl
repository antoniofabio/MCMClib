load ./mcmclib.so
source ./util.tcl

set lst [mcmclib_vector_list_alloc]
stopifne [vector_list_str_next_get $lst] NULL
