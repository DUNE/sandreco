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

iname=$(basename ${ifile})
idir=$(dirname ${ifile})

ddir="/wd/dune-it/enurec/analysis/files/digits"
rdir=$(echo ${ddir} | sed "s/digits/reco/g")
adir=$(echo ${rdir} | sed "s/reco/analysis/g")

if [ ! -f "${ifile}" ]
then
	 echo -e "\e[5m\e[91mERROR\e[0m: input file ${ifile} does not exist."
   return 1
fi

if [ ! -d "${ddir}" ]
then
	 echo -e "\e[5m\e[91mERROR\e[0m: digit directory ${ddir} does not exist. Please create it or change path in the script."
   return 1
fi

if [ ! -d "${rdir}" ]
then
	 echo -e "\e[5m\e[91mERROR\e[0m: reco directory ${rdir} does not exist. Please create it or change path in the script."
   return 1
fi

if [ ! -d "${rdir}" ]
then
	 echo -e "\e[5m\e[91mERROR\e[0m: reco directory ${rdir} does not exist. Please create it or change path in the script."
   return 1
fi

if [ ! -d "${adir}" ]
then
	 echo -e "\e[5m\e[91mERROR\e[0m: analysis directory ${adir} does not exist. Please create it or change path in the script."
   return 1
fi

dname=$(echo ${iname} | sed "s/edep-sim/digits/g")
rname=$(echo ${dname} | sed "s/digits/reco/g")
aname=$(echo ${rname} | sed "s/reco/analysis/g")

dfile="${ddir}/${dname}"
rfile="${rdir}/${rname}"
afile="${adir}/${aname}"

echo "===============================" &&
echo "Digitization                  =" && 
echo "===============================" && 
Digitize "${ifile}" "${dfile}"

echo "===============================" &&
echo "Reconstruction                =" && 
echo "===============================" && 
Reconstruct "${dfile}"

echo "===============================" &&
echo "Analysis                      =" && 
echo "===============================" && 
Analyze "${dfile}"