#! /bin/sh

# Always edit the project.cmt file for nightly, since the release might change
echo "Making sure your project.cmt file is up to date..."
cd cmt
rm -f project.cmt
echo "project MINERVA_$tag" >> project.cmt
echo "use MINERVA MINERVA_$MINERVA_RELEASE" >> project.cmt
setenvMinerva $tag

