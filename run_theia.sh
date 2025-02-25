#! /usr/bin/env bash

##################################################

HORN=$1
SWAP=$2
NPER=$3
STT=$4
if [ "${HORN}" != "FHC" ] && [ "${HORN}" != "RHC" ]; then
echo "Invalid beam mode ${HORN}"
echo "Must be FHC or RHC"
kill -INT $$
fi

MODE="neutrino"
MODE="${MODE}"
if [ "${HORN}" = "RHC" ]; then
MODE="anti${MODE}"
fi

CP="ifdh cp"
echo "Running edepsim for ${HORN} mode, ${NPER} events"

RNDSEED=$((${PROCESS}))

USERDIR="/pnfs/dune/persistent/users/gyang/end"
OUTDIR="/pnfs/dune/scratch/users/gyang/end"

##################################################

## Setup UPS and required products

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup dk2nu        v01_05_01b   -q e15:prof
setup genie        v2_12_10c    -q e15:prof
setup genie_xsec   v2_12_10     -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10     -q dkcharmtau
setup root v6_26_06a -q e20:p3913:prof
setup geant4 v4_11_0_p01c -q debug:e20:qt
setup ifdhc

setup cmake v3_22_0
export CXX=`which g++` # this might be specific for Fermilab?
export CC=`which gcc` # this might be specific for Fermilab?

##################################################

echo "${CP} ${USERDIR}/ratpac-two.tar.gz ${PWD}/ratpat-two.tar.gz"
${CP} ${USERDIR}/ratpac-two.tar.gz ${PWD}/ratpac-two.tar.gz
tar -xzf ratpac-two.tar.gz
cd ratpac-two
ls -ltr
mkdir build
cd build
rm -rf *
cmake -DCMAKE_INSTALL_PREFIX=. ..
make
cd ..
cd ..
. ratpac-two/ratpac.sh

${CP} ${USERDIR}/EndSimulations.tar.gz ${PWD}/EndSimulations.tar.gz
tar -xzf EndSimulations.tar.gz
cd EndSimulations
ls -ltr
mkdir build
cd build
rm -rf *
cmake -DCMAKE_INSTALL_PREFIX=. ..
make
cd ..
cd ..
. EndSimulations/end.sh

echo "for i in {$((${STT}))..$((STT+NPER))}"
for i in {$((${STT}))..$((STT+NPER))}
do
  echo "in the loop"
  echo "ifdh ls ${USERDIR}/genie_${MODE}_cc2023_${SWAP}_${STT}_$((${STT}+1)).root into ls1.txt"
  ifdh ls ${USERDIR}/genie_${MODE}_cc2023_${SWAP}_${STT}_$((${STT}+1)).root > ls1.txt
  if [ -s ls1.txt ]; then
    echo "${CP} ${USERDIR}/genie_${MODE}_cc2023_${SWAP}_${STT}_$((${STT}+1)).root EndSimulations/genie_${MODE}_cc2023_${SWAP}_current.root"
    ${CP} ${USERDIR}/genie_${MODE}_cc2023_${SWAP}_${STT}_$((${STT}+1)).root EndSimulations/genie_${MODE}_cc2023_${SWAP}_current.root
  else
    echo "file not found: ${USERDIR}/genie_${MODE}_cc2023_${SWAP}_${STT}_$((${STT}+1)).root"
  fi
  echo "end end.mac -o test1.root"
  cd EndSimulations

  cp end.mac end_temp.mac
  sed -i "s/neutrino/${MODE}/g" end_temp.mac
  sed -i "s/swap1/${SWAP}/g" end_temp.mac
  end end_temp.mac -o test1.root
  rm end_temp.mac

  echo "${CP} test1.root ${OUTDIR}/test1_${STT}_$((${STT}+1)).root"
  ${CP} test1.root ${OUTDIR}/test1_${STT}_$((${STT}+1)).root
  cd ..
done
##################################################


