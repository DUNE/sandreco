#!/bin/bash

source $(dirname ${BASH_SOURCE[0]})/env.sh

apps=(Digitize Reconstruct Analyze)

for a in ${apps[*]}
do 
   res=$(which ${a} 2> /dev/null)

   if [ "${res}" == "" ] 
   then
      echo -e "\e[5m\e[91mERROR\e[0m: ${a} not exist. Please set the right environment or build the executables."
      return 1
   fi
done

ifile=${1}

if [ "${ifile}" == "" ]
then
	echo -e "\e[5m\e[91mERROR\e[0m: source ${BASH_SOURCE[0]} <input file>"
	return 1
fi

iname=$(basename ${ifile})
idir=$(dirname ${ifile})

rdir="/home/dune-it/data/reco"

if [ ! -f "${ifile}" ]
then
	 echo -e "\e[5m\e[91mERROR\e[0m: input file ${ifile} does not exist."
   return 1
fi

if [ ! -d "${rdir}" ]
then
	 echo -e "\e[5m\e[91mERROR\e[0m: digit directory ${rdir} does not exist. Please create it or change path in the script."
   return 1
fi

rname=$(echo ${iname} | sed "s/edep-sim/reco/g")

rfile="${rdir}/${rname}"

echo "===============================" &&
echo "Digitization                  =" && 
echo "===============================" && 
Digitize "${ifile}" "${rfile}"

echo "===============================" &&
echo "Reconstruction                =" && 
echo "===============================" && 
Reconstruct "${rfile}"

echo "===============================" &&
echo "Analysis                      =" && 
echo "===============================" && 
Analyze "${rfile}"
