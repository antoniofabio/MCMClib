source("rawgsl.R")

mcmclib <- new.env(parent=rawgsl)
dyn.import(pst("mcmclib_",
               c(pst("mh_", c("alloc", "free", "update", "update_N")),
                 "gauss_rw_alloc", "gauss_mrw_alloc",
                 pst("amh_", c("alloc", "free", "update", "update_N")),
                 "gauss_scalar_am_alloc", "gauss_am_alloc")),
           mcmclib)
