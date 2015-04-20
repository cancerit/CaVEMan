#!/bin/bash

##########LICENCE##########
# Copyright (c) 2014-2015 Genome Research Ltd.
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

SOURCE_HTSLIB="https://github.com/samtools/htslib/archive/1.2.tar.gz"

done_message () {
    if [ $? -eq 0 ]; then
        echo " done."
        if [ "x$1" != "x" ]; then
            echo $1
        fi
    else
        echo " failed.  See setup.log file for error messages." $2
        echo "    Please check INSTALL file for items that should be installed by a package manager"
        exit 1
    fi
}

get_distro () {
  EXT=""
  DECOMP="gunzip -f"
  echo "$1"
  if [[ $2 == *.tar.bz2* ]] ; then
    EXT="tar.bz2"
    DECOMP="bzip2 -fd"
  elif [[ $2 == *.tar.gz* ]] ; then
    EXT="tar.gz"
  else
    echo "I don't understand the file type for $1"
    exit 1
  fi
  if hash curl 2>/dev/null; then
    curl -sS -o $1.$EXT -L $2
  else
    wget -nv -O $1.$EXT $2
  fi
  mkdir -p $1
  `$DECOMP $1.$EXT`
  tar --strip-components 1 -C $1 -xf $1.tar
}

if [ "$#" -ne "1" ] ; then
  echo "Please provide an installation path  such as /opt/pancan"
  exit 0
fi

INST_PATH=$1

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

# re-initialise log file
echo > $INIT_DIR/setup.log

# log information about this system
(
    echo '============== System information ===='
    set -x
    lsb_release -a
    uname -a
    sw_vers
    system_profiler
    grep MemTotal /proc/meminfo
    set +x
    echo; echo
) >>$INIT_DIR/setup.log 2>&1

# cleanup inst_path
mkdir -p $INST_PATH/bin
cd $INST_PATH
INST_PATH=`pwd`
cd $INIT_DIR

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

cd $SETUP_DIR

echo -n "Building htslib ..."
if [ -e $SETUP_DIR/htslib.success ]; then
  echo -n " previously installed ...";
else
  cd $SETUP_DIR
  (
  set -xe
  if [ ! -e htslib ]; then
    get_distro "htslib" $SOURCE_HTSLIB
  fi
  patch htslib/cram/cram_index.c < $INIT_DIR/patches/htslibcramindex.diff
  make -C htslib -j$CPU
  touch $SETUP_DIR/htslib.success
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build htslib."

export HTSLIB="$SETUP_DIR/htslib"

echo -n "Building CaVEMan ..."
if [ -e "$SETUP_DIR/caveman.success" ]; then
  echo -n " previously installed ...";
else
  cd $INIT_DIR
  (
    set -xe
    mkdir -p $INIT_DIR/c/bin
    make -j$CPU
    cp $INIT_DIR/bin/caveman $INST_PATH/bin/.
    make clean
    touch $SETUP_DIR/caveman.success
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build CaVEMan."

#add bin path for install tests
export PATH="$INST_PATH/bin:$PATH"

# cleanup all junk
rm -rf $SETUP_DIR

echo
echo
