#!/bin/bash

VERSION=$(cat version)

#make version header file
echo "#define PACKAGE_VERSION $VERSION" > src/version.h

#Update documentation
rm -r html/doc
/Applications/Doxygen.app/Contents/Resources/doxygen Doxyfile
# Add logo to manual
awk 'NR==43{print "\\includegraphics{../../img/logo.pdf}\\\\"}1' html/doc/latex/refman.tex > html/doc/latex/temp.tex
awk 'NR==45{print "\\url{http://bio.math.berkeley.edu/eXpress}\\\\"}1' html/doc/latex/temp.tex > html/doc/latex/refman.tex
make -C html/doc/latex
mv html/doc/html/* html/doc
rmdir html/doc/html
mv html/doc/latex/refman.pdf html/doc/express-doc.pdf
rm -r html/doc/latex


#Create new download directories
mkdir html/downloads/express-$VERSION
mkdir html/downloads/express-$VERSION/express-$VERSION-macosx
mkdir html/downloads/express-$VERSION/express-$VERSION-linux
mkdir html/downloads/express-$VERSION/express-$VERSION-src

#Populate download directories
#macosx
cp osx_build/src/express html/downloads/express-$VERSION/express-$VERSION-macosx
cp README html/downloads/express-$VERSION/express-$VERSION-macosx/README
cp LICENSE html/downloads/express-$VERSION/express-$VERSION-macosx/LICENSE
#linux
cp linux_build/src/express html/downloads/express-$VERSION/express-$VERSION-linux
cp README html/downloads/express-$VERSION/express-$VERSION-linux/README
cp LICENSE html/downloads/express-$VERSION/express-$VERSION-linux/LICENSE
#src
cp -R src html/downloads/express-$VERSION/express-$VERSION-src/src
cp README html/downloads/express-$VERSION/express-$VERSION-src/README
cp LICENSE html/downloads/express-$VERSION/express-$VERSION-src/LICENSE

#Tar download directories
cd html/downloads/express-$VERSION
tar -czf express-$VERSION-macosx.tgz express-$VERSION-macosx 
tar -czf express-$VERSION-linux.tgz express-$VERSION-linux
tar -czf express-$VERSION-src.tgz express-$VERSION-src
rm -r express-$VERSION-macosx
rm -r express-$VERSION-linux
rm -r express-$VERSION-src

#add new files to git
git add .

#return to start directory
cd ../../../


