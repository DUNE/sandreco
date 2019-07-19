#!/bin/bash

DIR=$(dirname ${BASH_SOURCE[0]})

if [[ -z "${ROOT_INCLUDE_PATH}" ]]; then
  export ROOT_INCLUDE_PATH=${DIR}/../struct
else
  export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH}:${DIR}/../struct
fi

export PATH=${PATH}:${DIR}/../bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${DIR}/../struct
