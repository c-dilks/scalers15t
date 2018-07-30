#!/bin/bash

ls datfiles/run*.dat | sed 's/^.*run//;s/_.*$//' > listaHave
awk '{print $1}' good_trans.dat > listaWant
#vimdiff lista*

# build list of runs which we need to download
> list.toDownload.txt
while read run; do
  if [ -z `grep $run listaHave` ]; then
    echo $run >> list.toDownload.txt
  fi
done < listaWant

# build list of runs which we have datfiles for, but don't want
# to be included in the analysis (you should move them to some
# subdirectory of datfiles, or remove them)
> list.toOmit.txt
while read run; do
  if [ -z `grep $run listaWant` ]; then
    echo $run >> list.toOmit.txt
  fi
done < listaHave

# remove temp files
rm lista{Have,Want}
