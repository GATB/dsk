#!/bin/bash
#--------------------------------------------------------------#
#         Continuous integration script for Jenkins            #
#--------------------------------------------------------------#
#
# Default mode :
# This script will exit with error (exit code 1) if any of its steps fails.
# To change this behaviour, choose DO_NOT_STOP_AT_ERROR in Jenkins (see below).
#--------------------------------------------------------------#
set +xv

echo "
-----------------------------------------
 Miscellaneous information 
-----------------------------------------
date      : `date`
hostname  : `hostname`
pwd       : `pwd`

-----------------------------------------
 Jenkins build parameters (user defined)
-----------------------------------------
BRANCH_TO_BUILD      : ${BRANCH_TO_BUILD}
INRIA_FORGE_LOGIN    : ${INRIA_FORGE_LOGIN}
DO_NOT_STOP_AT_ERROR : ${DO_NOT_STOP_AT_ERROR}

-----------------------------------------
 Jenkins build parameters (built in)
-----------------------------------------
BUILD_NUMBER         : ${BUILD_NUMBER}
"

error_code () { [ "$DO_NOT_STOP_AT_ERROR" = "true" ] && { return 0 ; } }

[ "$DO_NOT_STOP_AT_ERROR" != "true" ] && { set -e ; } || { echo "(!) DEBUG mode, the script will NOT stop..." ; echo; }
set -xv

# quick look at resources
#-----------------------------------------------
free -h
#-----------------------------------------------
lstopo
#-----------------------------------------------
df -kh
#-----------------------------------------------


################################################################
#                       COMPILATION                            #
################################################################

gcc --version
g++ --version

[ `gcc -dumpversion` = 4.7 ] && { echo "GCC 4.7"; } || { echo "GCC version is not 4.7, we exit"; exit 1; }

JENKINS_TASK=tool-${TOOL_NAME}-build-debian7-64bits-gcc-4.7
GIT_DIR=/scratchdir/builds/workspace/gatb-${TOOL_NAME}
BUILD_DIR=/scratchdir/$JENKINS_TASK/gatb-${TOOL_NAME}/build

rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR

#-----------------------------------------------
# we need gatb-core submodule to be initialized
cd $GIT_DIR
git submodule init
git submodule update

#-----------------------------------------------
cd $BUILD_DIR

#-----------------------------------------------
cmake -Wno-dev -DJENKINS_TAG=${BRANCH_TO_BUILD} $GIT_DIR

#-----------------------------------------------
make -j 2 || error_code

################################################################
#                       TEST                                   #
################################################################
# prepare data and scripts
cp -R $GIT_DIR/scripts/ ..
cp -R $GIT_DIR/test/ ..
# run tests
cd ../scripts
./simple_test.sh || error_code
# cleanup disk space
cd ..
rm -rf scripts test
# go bask to build for packaging step
cd build

################################################################
#                       PACKAGING                              #
################################################################

# Upload bin bundle to the forge
if [ $? -eq 0 ] && [ "$INRIA_FORGE_LOGIN" != none ] && [ "$DO_NOT_STOP_AT_ERROR" != true ]; then
	make package
    scp ${TOOL_NAME}-${BRANCH_TO_BUILD}-bin-Linux.tar.gz ${INRIA_FORGE_LOGIN}@scm.gforge.inria.fr:/home/groups/gatb-tools/htdocs/ci-inria
    # source package is handled by the osx task
fi

