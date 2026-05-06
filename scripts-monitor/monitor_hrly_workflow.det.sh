#!/bin/bash
 
# current script path
script_dir="$(cd "$(dirname "$0")" && pwd)"

# load the rocoto module
source /etc/profile.d/modules.sh
ROCOTOMODULE=/scratch4/BMC/zrtrr/gge/rocoto/modulefiles
module use ${ROCOTOMODULE}
module load rocoto/1.3.7g

# Go to the real-time run
cd /scratch3/BMC/rtwrfruc/RRFSv2X/NCO_dirs/RRFSv2X/exp/rrfsdet


recipients="Liaofan.Lin@noaa.gov Liaofan.Lin@noaa.gov"
subject="rrfsv2x alert - det "

curtime=$(date -u +%Y%m%d%H)
PDY=${curtime:0:8}
cyc=${curtime:8:2}

#prevCyc=$(date -u -d "${PDY} ${cyc} UTC -3 hours" +%Y%m%d%H) # check the recent 3 cycles
prevCyc=$(date -u -d "${PDY} ${cyc} UTC -6 hours" +%Y%m%d%H) # check the recent 6 cycles
cycles=${prevCyc}00:${PDY}${cyc}00

send_email=false
msg=$(rocotostat -w rrfs.xml -d rrfs.db -c ${cycles})
echo "${msg}" > ${script_dir}/detmsg.new

# need to send out alerts when we get DEAD jobs
if [[ -s ${script_dir}/detmsg.save  && -s ${script_dir}/detmsg.new ]]; then
  # only send alerts if the new msg is different from the saved one to avoid duplicate alerts
  if ! diff ${script_dir}/detmsg.save ${script_dir}/detmsg.new &>/dev/null && [[ ${msg} == *DEAD* ]]; then
    send_email=true
  fi
elif [[ ${msg} == *DEAD* ]]; then
  send_email=true
fi
mv ${script_dir}/detmsg.new ${script_dir}/detmsg.save

## check whether if the workflow is stalled, i.e. no running jobs
## and if stalled for more than an hour, send out alerts
if [[ ${msg} != *RUNNING* ]]; then
  if [[ -s ${script_dir}/detstall ]]; then
    stall_bgn=$(head -n 1 ${script_dir}/detstall)
    cur_seconds=$(date +%s)
    diff_secons=$(( cur_seconds - stall_bgn ))
    if (( diff_secons > 3600 )); then # stalled for an hour
      send_email=true
      rm -rf ${script_dir}/detstall
    fi
  else
    date +%s > ${script_dir}/detstall
  fi
fi

if ${send_email}; then
  echo "${msg}" | mail -s "${subject}" ${recipients}
fi
