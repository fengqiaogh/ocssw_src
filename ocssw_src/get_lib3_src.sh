#! /bin/bash

#
# This script is used to grab the lib3 sources.
#
# $OCSSWROOT/OCSSW_bash.env must be sourced first
#

LIB3_TAG=`cat $OCSSWROOT/.manifest_tag`

$OCSSWROOT/src/manifest/install_ocssw.py --tag $LIB3_TAG --opt_src
