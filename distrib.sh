make clean
make doc
VNAME=`git describe --tags --abbrev=5`
git archive --format=tar --prefix=MCMClib_$VNAME/ HEAD | gzip >MCMClib_$VNAME.tar.gz
