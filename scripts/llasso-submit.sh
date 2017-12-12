#!/bin/sh

# Initial submit file for llasso-script.jl

module load R
module load julia
export LD_LIBRARY_PATH=`R RHOME`/lib

PRJDIR="${HOME}/22q"
#DATADIR="${PRJDIR}/data"
#OUTDIR="${PRJDIR}/output/FastQC"

if [ -e /bin/mktemp ]; then
	TMPDIR=`/bin/mktemp -d /scratch/XXXXXX`
elif [ -e /usr/bin/mktemp ]; then
	TMPDIR=`/usr/bin/mktemp -d /scratch/XXXXXX`
else
	echo "Error. Cannot find program to create tmp directory"
	exit
fi

cp ${PRJDIR}/* ${TMPDIR}
cd ${TMPDIR}

julia ${TMPDIR}/llasso-script.jl

/bin/rm ${TMPDIR}/*.fam
/bin/rm ${TMPDIR}/*.bim
/bin/rm ${TMPDIR}/*.bed

rsync -av ${TMPDIR}/ ${PRJDIR}

/bin/rm -fr ${TMPDIR}

module unload R
module unload julia