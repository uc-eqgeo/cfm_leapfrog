#!/bin/sh

alias trelis=/Applications/Trelis-17.1.app/Contents/MacOS/Trelis-17.1
for jou in *.jou
do
  trelis -nographics -nojournal "$jou"
  stl_name="${jou%.jo*}_volume.stl"
  obj_name="${jou%.jo*}_volume.obj"
  meshio convert "${stl_name}" "${obj_name}"
done

rm volumes.zip
zip volumes.zip *.obj