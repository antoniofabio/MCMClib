make clean
VNAME=`date '+%m-%d-%Y'`
git archive --format=tar --prefix=MCMClib_$VNAME/ HEAD | gzip >MCMClib_$VNAME.tar.gz
