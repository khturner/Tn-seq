#!/bin/bash

usage () {
  echo "usage: $0 [-i <#>] [-e <#>] [-a <assembly>] [-o <output>] [-c <name>] [-x <#>] <pfx1> <pfx2> <pfx3> ... <pfxn> "
  echo "Required parameters:"
  echo "-i     The number of the most represented sites to ignore"
  echo "-e     The number of randomly distributed pseudo-datasets to generate"
  echo "-a     The name of the assembly you're using (PAO1, PA14, AAD7S, SGCH1, ECK12W3110, HG003)"
  echo "-o     The name for the output file"
  echo "-c     The name for the control condition"
  echo "-x     The number of replicates for the control condition"
  echo ""
  echo "The required parameters must precede the prefixes to be joined, listed with the"
  echo "  control conditions followed by the test conditions. See examples below."
  echo ""
  echo "Examples:"
  echo "$0 -i 50 -e 10 -a PAO1 -o Example -c control -x 2 C1 C2"
}

# Read in the important options
while getopts ":i:e:a:o:h:c:x:" option; do
  case "$option" in
  	i)  CUT="$OPTARG" ;;
  	e)	PSEUDO="$OPTARG" ;;
  	a)  ASSEMBLY="$OPTARG" ;;
    o)  OUT_PFX="$OPTARG" ;;
    c)  CONTROL_PFX="$OPTARG" ;;
    x)  CONTROL_REPS="$OPTARG" ;;
    h)  # it's always useful to provide some help 
        usage
        exit 0 
        ;;
    :)  echo "Error: -$option requires an argument" 
        usage
        exit 1
        ;;
    ?)  echo "Error: unknown option -$option" 
        usage
        exit 1
        ;;
  esac
done    
shift $(( OPTIND - 1 ))

# Do some error checking to make sure parameters are defined
if [ -z "$CUT" ]; then
  echo "Error: you must specify the number of sites to ignore using -i"
  usage
  exit 1
fi

if [ -z "$PSEUDO" ]; then
  echo "Error: you must specify the number of pseudo-datasets to generate using -e"
  usage
  exit 1
fi

if [ -z "$ASSEMBLY" ]; then
  echo "Error: you must specify an assembly using -a"
  usage
  exit 1
fi

if [ -z "$OUT_PFX" ]; then
  echo "Error: you must specify an output prefix for your file using -o"
  usage
  exit 1
fi

if [ -z "$CONTROL_PFX" ]; then
  echo "Error: you must specify the name for your control"
  echo "using -c"
  usage
  exit 1
fi

if [ -z "$CONTROL_REPS" ]; then
  echo "Error: you must specify the number of replicates for your control"
  echo "condition using -x"
  usage
  exit 1
fi

# Give the usage if there aren't enough parameters
if [ $# -lt 1 ] ; then
  echo "Please provide at least 1 file"
  usage
  exit 1
fi

ASSEMBLY_PFX="/home1/02127/khturner/ref_genome/$ASSEMBLY/$ASSEMBLY"

# Smoothing (LOESS) and normalization (TMM)
echo "Performing LOESS smoothing, normalization and obligate essentiality analysis on count data..."
R --vanilla --args $CONTROL_PFX $CONTROL_REPS $ASSEMBLY_PFX $OUT_PFX $CUT $PSEUDO $@ < ~/local/bin/TnSeqDESeqEssential.R