#!/bin/bash

root="root"

res=$(which ${root} 2> /dev/null)

if [ "${res}" == "" ] 
then
	 echo -e "\e[5m\e[91mERROR\e[0m: ${root} not exist. Please set the right environment."
   return 1
fi

dlib="/wd/dune-it/enurec/analysis/kloe-simu/digitization/digitization_cpp.so"
rlib="/wd/dune-it/enurec/analysis/kloe-simu/reconstruction/reconstruction_cpp.so"
alib="/wd/dune-it/enurec/analysis/kloe-simu/analysis/analysis_cpp.so"

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

if [ ! -f "${dlib}" ]
then
	 echo -e "\e[5m\e[91mERROR\e[0m: ${dlib} does not exist. Please create it or change path in the script."
   return 1
fi

if [ ! -f "${rlib}" ]
then
	 echo -e "\e[5m\e[91mERROR\e[0m: ${rlib} does not exist. Please create it or change path in the script."
   return 1
fi

if [ ! -f "${alib}" ]
then
	 echo -e "\e[5m\e[91mERROR\e[0m: ${alib} does not exist. Please create it or change path in the script."
   return 1
fi

dname=$(echo ${iname} | sed "s/edep-sim/digits/g")
rname=$(echo ${dname} | sed "s/digits/reco/g")
aname=$(echo ${rname} | sed "s/reco/analysis/g")

dfile="${ddir}/${dname}"
rfile="${rdir}/${rname}"
afile="${adir}/${aname}"

echo "gSystem->Load(\"${dlib}\"); Digitize(\"${ifile}\",\"${dfile}\");" | root -l
echo "gSystem->Load(\"${rlib}\"); Reconstruct(\"${dfile}\",\"${ifile}\",\"${rfile}\");" | root -l
echo "gSystem->Load(\"${alib}\"); Analyze(\"${rfile}\",\"${ifile}\",\"${afile}\");" | root -l

#echo "gSystem->Load(\"/wd/dune-it/enurec/analysis/kloe-simu/reconstruction/reconstruction_cpp.so\"); Reconstruct(\"/wd/dune-it/enurec/analysis/files/digits/numu_geoV12_1000.0.digits.root\",\"/wd/dune-it/enurec/prod/data/edep-sim/numu_geoV12_1000.0.edep-sim.root\",\"/wd/dune-it/enurec/analysis/files/reco/numu_geoV12_1000.0.reco.root\")" | root -l

#echo "gSystem->Load(\"/wd/dune-it/enurec/analysis/kloe-simu/analysis/analysis_cpp.so\"); Analyze(\"/wd/dune-it/enurec/analysis/files/reco/numu_geoV12_1000.0.reco.root\",\"/wd/dune-it/enurec/prod/data/edep-sim/numu_geoV12_1000.0.edep-sim.root\",\"/wd/dune-it/enurec/analysis/files/analysis/numu_geoV12_1000.0.analysis.root\")" | root -l
