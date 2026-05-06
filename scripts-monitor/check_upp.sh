#!/usr/bin/env bash

source /etc/profile

script_dir="$(cd "$(dirname "$0")" && pwd)"


cd /scratch3/BMC/rtwrfruc/RRFSv2X/NCO_dirs/RRFSv2X/exp/rrfsdet
source  qrocoto/load_qrocoto.sh


#  module load rocoto
#rtlog="/scratch4/BMC/rtrr/NCO_dirs/code/rtlog.rrfsv2x"

#while true; do

  rtasks upp_g01 36  > ${script_dir}/rtlog-upp.rrfsv2x  &&  mail -s 'Ursa: RRFSv2X status upp check' liaofan.lin@noaa.gov < ${script_dir}/rtlog-upp.rrfsv2x
  
  #  sleep 43200

#done
