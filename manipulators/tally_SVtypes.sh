#!/bin/bash

# Simple function to tally up and pretty print (to stdout) SV types from SMAP (v0.7) 

if [ $# -ne 1 ]
then
	echo "usage: tally_SVtypes input.smap"
	exit
fi
if [ ! -e $1 ] 
then
	echo "Can't find ${1}"
	exit 
fi

# Assumes Type is the 10th column of input smap file!!
mytypes=$(grep -v "^#" ${1} | cut -f10)
array=(${mytypes})

# Sum the types
ndel=0; ndup=0; nins=0; ninv=0; nintra=0; ninter=0; ninvp=0;
for type in ${array[*]}; do
	case $type in 
		deletion*)
			(( ndel++ ))
			;;
		duplication*)
			(( ndup++ ))
			;;
		insertion*)
			(( nins++ ))
			;;
		inversion_partial)
			(( ninv++ ))
			;;
		inversion_paired)
			(( ninvp++ ))
			;;
		*intrachr*)
			(( nintra++ ))
			;;
		*interchr*)
			(( ninter++ ))
			;;
	esac
done
# inversions are paired: each inversion and inversion_nbase has a partner of type inversion_partial; while each paired-inversion will have two smap entries both of type inversion_paired
ninv=$(( $ninv + (( $ninvp / 2 )) ))
# total: 
ntot=$(( $nins + $ndel + $ndup + $ninv + $nintra + $ninter ))

echo -e "insertions\t${nins}"
echo -e "deletions\t${ndel}"
echo -e "duplications\t${ndup}"
echo -e "inversions\t${ninv}"
echo -e "intra-chr\t${nintra}"
echo -e "inter-chr\t${ninter}"
echo -e "TOTAL SV\t$ntot"




