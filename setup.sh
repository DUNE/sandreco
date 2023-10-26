#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# SANDRECO
export PATH=${DIR}/bin:${PATH}
export LD_LIBRARY_PATH=${DIR}/lib:${LD_LIBRARY_PATH}
export sandreco_DIR=${DIR}
