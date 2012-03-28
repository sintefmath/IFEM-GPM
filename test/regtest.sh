#!/bin/bash

mysim=$1

cd `dirname $2`
readarray < $2
$mysim $MAPFILE 2>&1 > templog
globres=1
IFS=$'\n'

resultfile=`basename $2 | sed 's/\.reg/\.res/'`
echo $resultfile

difference=`diff templog $resultfile`
if ! test $? -eq 0
then
  globres=0
elif test -z $difference
then
  globres=1
else
  globres=0
fi

echo "globres=$globres"

if test $globres -eq 0
then
  cat templog >> ../../../failed.log
  rm templog
  exit 1
fi

rm templog
exit 0
