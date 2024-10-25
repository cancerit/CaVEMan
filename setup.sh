#!/bin/bash

##########LICENCE##########
# Copyright (c) 2014-2016 Genome Research Ltd.
#
# Author: CancerIT <cgpit@sanger.ac.uk>
#
# This file is part of CaVEMan.
#
# CaVEMan is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
##########LICENCE##########

REQUIRED_MIN_LIBZ="1.2.3.4"

function version_eq_gt() {
    if [ "$1" = "$2" ]; then
        return 0
    fi
    test "$(printf '%s\n' "$@" | sort -V | head -n 1)" != "$1";
}

set -euo pipefail 

if [ "$#" -lt "1" ] ; then
  echo "Please provide an installation path  such as /opt/pancan"
  exit 0
fi

INST_PATH=$1
DEBUG_ARG="${2:-}"
if [ "$DEBUG_ARG" = "DEBUG" ]; then
  DEBUG_FLAG=1
  echo "DEBUG set"
else
  DEBUG_FLAG=0
fi

CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi
echo "Max compilation CPUs set to $CPU"

# get current directory
INIT_DIR=`pwd`

echo '#include <zlib.h>' | cpp -H 1>/dev/null 2> tmp.log.cpp
LIBZ_VER=$(cat tmp.log.cpp | head -1 | cut -d' ' -f 2 | xargs grep -e '#define ZLIB_VERSION' | cut -d ' ' -f 3 | sed 's/"//g')
rm tmp.log.cpp

if version_eq_gt $LIBZ_VER $REQUIRED_MIN_LIBZ ; then
	echo "Found acceptable libz version $LIBZ_VER."
	echo "Continuing install"
else
	echo "ERROR: CaVEMan requires libz version >= $REQUIRED_MIN_LIBZ"
	echo "Found libz version: $LIBZ_VER"
	echo "Exiting install"
	exit 1
fi

# cleanup inst_path
mkdir -p $INST_PATH/bin
cd $INST_PATH
INST_PATH=`pwd`
cd $INIT_DIR

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

cd $SETUP_DIR
echo -n "Building linasm..."
if [ -e $SETUP_DIR/linasm.success ]; then
  echo -n "previously installed..."
else
  wget https://github.com/rurban/linasm/archive/e33b4a083f8bbdcbd58018c9abc65047d20fd431.zip &&
  unzip -o e33b4a083f8bbdcbd58018c9abc65047d20fd431.zip && mv linasm-e33b4a083f8bbdcbd58018c9abc65047d20fd431 linasm && cd linasm &&
  mkdir -p linasm_inst &&
  make &&
  make install prefix=${SETUP_DIR}/linasm/linasm_inst/ &&
  make clean &&
  touch $SETUP_DIR/linasm.success
fi

export LINASM_INC=$SETUP_DIR/linasm/linasm_inst/include/
export LINASM_LIB=$SETUP_DIR/linasm/linasm_inst/lib/

cd $SETUP_DIR
echo -n "Building htslib ..."
if [ -e $SETUP_DIR/htslib.success ]; then
  echo -n "previously installed ...";
else
  wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 &&
  bzip2 -fd htslib-1.10.2.tar.bz2 &&
  mkdir -p htslib &&
  tar --strip-components 1 -C htslib -xf htslib-1.10.2.tar &&
  cd htslib &&
  make -j$CPU &&
  touch $SETUP_DIR/htslib.success
fi

export HTSLIB="$SETUP_DIR/htslib"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${HTSLIB}:${LINASM_LIB}"

echo -n "Building CaVEMan ..."
cd $INIT_DIR
make clean &&
{
if [ $DEBUG_FLAG -eq 1 ]; then
  make -j$CPU DEBUG=1
else
  make -j$CPU
fi
} &&
mv $INIT_DIR/bin/* $INST_PATH/bin/ &&
mkdir -p $INST_PATH/lib/ &&
cp $LINASM_LIB/liblinasm.so $INST_PATH/lib/. &&
touch $SETUP_DIR/caveman.success

# cleanup all intermediates
if [ -e "$SETUP_DIR/caveman.success" ]; then
  echo
  echo "Installation succesful"
  echo "Binaries available at $INST_PATH/bin/"
  echo "linasm.so available at $INST_PATH/lib/ - directory must be on LD_LIBRARY_PATH"
  if [ $DEBUG_FLAG -eq 1 ]; then
    echo "DEBUG set, retaining intermediate files"
  else
    rm -rf $SETUP_DIR
    make clean
  fi
else
  echo
  echo "Installation failed"
  echo "Retaining intermediate files"
fi

