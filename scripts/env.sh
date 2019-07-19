#!/bin/bash

if [[ -z "${ROOT_INCLUDE_PATH}" ]]; then
  export ROOT_INCLUDE_PATH=${BASH_SOURCE[0]}/../struct
else
  export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH}:${BASH_SOURCE[0]}/../struct
fi

export PATH=${PATH}:${BASH_SOURCE[0]}/../bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${BASH_SOURCE[0]}/../struct
