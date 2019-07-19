#!/bin/bash

if [[ -z "${ROOT_INCLUDE_PATH}" ]]; then
  export ROOT_INCLUDE_PATH=/wd/dune-it/enurec/analysis/kloe-simu/struct
else
  export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH}:/wd/dune-it/enurec/analysis/kloe-simu/struct
fi

export PATH=${PATH}:/wd/dune-it/enurec/analysis/kloe-simu/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/wd/dune-it/enurec/analysis/kloe-simu/struct
