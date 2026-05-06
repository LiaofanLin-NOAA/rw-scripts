#!/usr/bin/env bash

source /etc/profile

cd /scratch3/BMC/rtwrfruc/RRFSv2X/NCO_dirs/RRFSv2X/exp/v2.1.3/rrfsdet
source  qrocoto/load_qrocoto.sh
#  module load rocoto
rtlog="/scratch4/BMC/rtrr/NCO_dirs/code/rtlog.rrfsv2x"

#while true; do

  rtasks upp_g01 36  > ${rtlog}  &&  mail -s 'Ursa: RRFSv2X status check' chunhua.zhou@noaa.gov < ${rtlog}
#  sleep 43200

#done
