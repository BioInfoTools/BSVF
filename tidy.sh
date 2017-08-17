#!/bin/bash

PLFiles="bsuit lib/BSuitInc.pm lib/BSuitLib.pm"
COMMOMTIDY="-b -l=200 -t -nola"

perltidy $COMMOMTIDY -dbc -dsc -dp $PLFiles
perltidy $COMMOMTIDY -mbl=1 $PLFiles

