# Benchmark scripts for Speck annotation checker [![](https://images.microbadger.com/badges/image/satta/speck-benchmarks.svg)](https://microbadger.com/images/satta/speck-benchmarks "Get your own image badge on microbadger.com")

`gal_test.pl` -- GAL/Bioperl variant

  - Requires GAL (http://www.sequenceontology.org/software/GAL.html) 
    and Bioperl (from CPAN)

`gffutils_test.py` -- gffutils/Biopython variant

  - Requires gffutils (https://pythonhosted.org/gffutils) and Biopython
    (both installable via pip)
  - Python3 compatible

specs directory

  - contains use cases

## Usage

```
docker run -v /tmp:/mnt/out satta/speck-benchmarks /opt/scripts/run_performance_checks.sh /opt/genomes /mnt/out /opt/scripts

```
