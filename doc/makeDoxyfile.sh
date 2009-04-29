VERSION=`git describe --tags --abbrev=5`
sed "s/VERSION-NUMBER/${VERSION}/g" Doxyfile.in > Doxyfile
