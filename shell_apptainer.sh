#!/bin/bash
if [[ $# -eq 0 ]]; then
  exec /cvmfs/uboone.opensciencegrid.org/bin/shell_apptainer.sh -d
else
  exec /cvmfs/uboone.opensciencegrid.org/bin/shell_apptainer.sh "$@"
fi
