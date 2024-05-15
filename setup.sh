# Running SL7 inside a container
# The following command initiates an Apptainer (previously known as Singularity) shell session 
# with a Scientific Linux 7 (SL7) image from Fermilab's repository. The command binds several 
# directories to the container to ensure they are accessible from within the container.

# Launching the Apptainer shell with SL7
/cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer shell --shell=/bin/bash \
-B /cvmfs,/exp,/nashome,/opt,/run/user,/etc/hostname,/etc/hosts,/etc/krb5.conf --ipc --pid \
/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:latest

# Or /pnfs/dune
# /cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer shell --shell=/bin/bash \
# -B /cvmfs,/exp,/nashome,/pnfs/dune,/opt,/run/user,/etc/hostname,/etc/hosts,/etc/krb5.conf --ipc --pid \
# /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:latest

# If you are planning to set up UPS (Unix Product Support) software inside the SL7 container 
# you should also set the UPS_OVERRIDE environment variable before doing any setup:
export UPS_OVERRIDE="-H Linux64bit+3.10-2.17"

# Navigate to the directory where additional setup will be done
cd /exp/dune/app/users/mazam/production/CAFs/MiniRun5

# Source the DUNE setup script to configure the environment for DUNE software
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

# Set the keyboard layout to US
setxkbmap us

# Setup various required software packages using UPS
setup cmake v3_22_2
setup gcc v9_3_0
setup pycurl
setup ifdhc
setup geant4 v4_11_0_p01c -q e20:debug
setup dk2nugenie   v01_10_01k -q debug:e20
setup genie_xsec   v3_04_00 -q AR2320i00000:e1000:k250
setup genie_phyopt v3_04_00 -q dkcharmtau
setup jobsub_client
setup eigen v3_3_5
setup duneanaobj v03_02_01 -q e20:prof
setup hdf5 v1_10_5a -q e20
setup fhiclcpp v4_15_03 -q debug:e20

# Edep-sim needs to know where a certain GEANT4 .cmake file is located...
G4_cmake_file=`find ${GEANT4_FQ_DIR}/lib -name 'Geant4Config.cmake'`
export Geant4_DIR=`dirname $G4_cmake_file`

# Edep-sim needs to have the GEANT4 bin directory in the PATH
export PATH=$PATH:$GEANT4_FQ_DIR/bin

# Suppress ROOT include errors
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$GENIE_INC/GENIE

# nusystematics paths
# export NUSYST=${PWD}/nusystematics
# export LD_LIBRARY_PATH=${NUSYST}/build/Linux/lib:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=${NUSYST}/build/nusystematics/artless:$LD_LIBRARY_PATH
# export FHICL_FILE_PATH=${NUSYST}/nusystematics/fcl:$FHICL_FILE_PATH

# Add pyGeoEff to PYTHONPATH
export PYTHONPATH=${PYTHONPATH}:${PWD}/DUNE_ND_GeoEff/lib/

# duneanaobj needs to be in the library paths too
export LD_LIBRARY_PATH=${DUNEANAOBJ_LIB}:$LD_LIBRARY_PATH

# Finally, add our lib & bin directories to the paths
mydir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export LD_LIBRARY_PATH=$mydir/lib:$LD_LIBRARY_PATH
export PATH=$mydir/bin:$PATH
