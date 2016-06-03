#!/bin/bash

cd $1
mkdir -p $2
for f in `ls | cut -f 1 -d . | sort | uniq`; do
  if [ -e $f.noseq.gff3.gz ]; then
    gunzip $f.noseq.gff3.gz
    rm -f $f.run.o $f.run.e
    rm -f $2/$f.runtime;
    $3/measure_space_peak.sh $3/gal_test.pl \
      $f.noseq.gff3 $f.genome.fasta > $f.gal.runtime.tmp 2>&1
    gffutils-cli create $f.noseq.gff3
    $3/measure_space_peak.sh $3/gffutils_test.py \
      $f.noseq.gff3.db $f.genome.fasta > $f.gffutils.runtime.tmp 2>&1
    $3/measure_space_peak.sh gt speck \
      -specfile specs/smallspec.lua \
      -seqfile $f.genome.fasta -matchdescstart \
      -typecheck sofa \
      -colored no -output statsonly -v \
      $f.noseq.gff3 > $f.gt.runtime.tmp 2>&1
     echo -ne \"$f\t\" > $2/$f.runtime
     grep time: $f.gal.runtime.tmp | cut -f 2 -d: | xargs echo -n >> $2/$f.runtime
     echo -ne \"\t\" >> $2/$f.runtime
     grep VmHWM: $f.gal.runtime.tmp | cut -f 5 -d' ' | xargs echo -n >> $2/$f.runtime
     echo -ne \"\t\" >> $2/$f.runtime
     grep time: $f.gffutils.runtime.tmp | cut -f 2 -d: | xargs echo -n >> $2/$f.runtime
     echo -ne \"\t\" >> $2/$f.runtime
     grep VmHWM: $f.gffutils.runtime.tmp | cut -f 5 -d' ' | xargs echo -n >> $2/$f.runtime
     echo -ne \"\t\" >> $2/$f.runtime
     grep time: $f.gt.runtime.tmp | cut -f 2 -d: | xargs echo -n >> $2/$f.runtime
     echo -ne \"\t\" >> $2/$f.runtime
     grep VmHWM: $f.gt.runtime.tmp | cut -f 5 -d' '| xargs echo -n >> $2/$f.runtime
     echo -ne \"\t\" >> $2/$f.runtime
     head -n 1 $f.gt.runtime.tmp >> $2/$f.runtime
  fi
done
