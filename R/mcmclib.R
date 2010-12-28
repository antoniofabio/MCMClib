source("rawgsl.R")

mcmclib <- new.env(parent=rawgsl)
dyn.import(pst("mcmclib_",
               c(pst("mh_", c("alloc", "free", "update", "update_N")),
                 "gauss_rw_alloc", "gauss_mrw_alloc",
                 pst("amh_", c("alloc", "free", "update", "update_N")),
                 "gauss_scalar_am_alloc", "gauss_am_alloc",
                 pst("monitor_", c("alloc", "free", "update")),
                 pst("monitor_get_", c("means", "vars", "ar", "msjd")),
                 pst("monitor_ecdf_", c("alloc", "free", "update")),
                 pst("monitor_acf_", c("alloc", "free", "update", "get")),
                 "iact_from_acf")),
           mcmclib)

with(mcmclib, {
  monitor_ecdf.type <- structType(list(Fn = pointerType,
                                       X0 = pointerType,
                                       n = doubleType,
                                       workspace = pointerType))
  monitor_ecdf_get_Fn <- function(m) {
    return(gvec2vec(getStructField(m, "Fn", monitor_ecdf.type)))
  }
})
