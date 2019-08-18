#!/bin/bash

DIR=$(dirname ${BASH_SOURCE[0]})

if [[ -z "${ROOT_INCLUDE_PATH}" ]]; then
  export ROOT_INCLUDE_PATH=${DIR}/../lib
else
  export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH}:${DIR}/../struct:/usr/include/eigen3
fi

export PATH=${PATH}:${DIR}/../bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${DIR}/../lib
