#!/usr/bin/env bash

usage () {
  echo "usage: $0 [-i <#>] [-a <assembly>] [-o <output>] [-c <name>] [-x <#>] [-t <name>] [-y <#>] <pfx1> <pfx2> <pfx3> ... <pfxn> "
  echo "Required parameters:"
  echo "-i     The number of the most represented sites to ignore"
  echo "-a     The name of the assembly you're using (PAO1 only so far)"
  echo "-o     The name for the output file"
  echo "-c     The name for the control condition"
  echo "-x     The number of replicates for the control condition"
  echo "-t     The name for the test condition"
  echo "-y     The number of replicates for the test condition"
  echo ""
  echo "The required parameters must precede the prefixes to be joined, listed with the"
  echo "  control conditions followed by the test conditions. See example below."
  echo ""
  echo "Example:"
  echo "$0 -i 50 -a PAO1 -o Example -c control -x 2 -t test -y 2 C1 C2 T1 T2"
}

# Read in the important options
while getopts ":i:a:o:h:c:x:t:y:" option; do
  case "$option" in
  	i)  CUT="$OPTARG" ;;
  	a)  ASSEMBLY="$OPTARG" ;;
    o)  OUT_PFX="$OPTARG" ;;
    c)  CONTROL_PFX="$OPTARG" ;;
    x)  CONTROL_REPS="$OPTARG" ;;
    t)  TEST_PFX="$OPTARG" ;;
    y)  TEST_REPS="$OPTARG" ;;
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

if [ -z "$TEST_REPS" ]; then
  echo "Error: you must specify the number of replicates for your test"
  echo "condition using -y"
  usage
  exit 1
fi

if [ -z "$TEST_PFX" ]; then
  echo "Error: you must specify the number of replicates for your test"
  echo "condition using -t"
  usage
  exit 1
fi

# Give the usage if there aren't enough parameters
if [ $# -lt 2 ] ; then
  echo "you cannot compare less than 2 files"
  usage
  exit 1
fi

ASSEMBLY_PFX="$REFGENOME/$ASSEMBLY/$ASSEMBLY"

# Trim most represented insertions and prep sites.txt file for analysis
echo "Trimming $CUT most represented sites and preparing sites.txt for analysis..."
for PREFIX in "$@"
do
	tail -n +$CUT $PREFIX/$PREFIX-sites.txt | sort -k 2 > $PREFIX/$PREFIX-tailsites.txt
	PREFIXES[$INDEX]=$PREFIX
	let "INDEX += 1"
done
join -a1 -a2 -1 2 -2 2 -o 0 1.1 2.1 -e "0" ${PREFIXES[0]}/${PREFIXES[0]}-tailsites.txt ${PREFIXES[1]}/${PREFIXES[1]}-tailsites.txt > $OUT_PFX-sites.txt
INDEX=2
while [ $INDEX -lt ${#PREFIXES[@]} ]
do
	COL=3
	FIELDS="1.2 1."$COL" "
	while [ $COL -le $INDEX ]
	do
		let "COL += 1"
		FIELDS=$FIELDS"1."$COL" "
	done
	join -a1 -a2 -1 1 -2 2 -o 0 $FIELDS 2.1 -e "0" $OUT_PFX-sites.txt ${PREFIXES[$INDEX]}/${PREFIXES[$INDEX]}-tailsites.txt > $OUT_PFX-newsites.txt
	cp $OUT_PFX-newsites.txt $OUT_PFX-sites.txt
	let "INDEX += 1"
done
echo "Position "${PREFIXES[@]} > $OUT_PFX-sites.txt
cat $OUT_PFX-newsites.txt >> $OUT_PFX-sites.txt
rm $OUT_PFX-newsites.txt

# Smoothing (LOESS) and normalization (TMM)
echo "Performing LOESS smoothing, normalization and differential abundance analysis on count data..."
R --vanilla --args $CONTROL_PFX $CONTROL_REPS $TEST_PFX $TEST_REPS $ASSEMBLY_PFX $OUT_PFX < ~/local/bin/TnSeqDESeq3.R > /dev/null 2> /dev/null
