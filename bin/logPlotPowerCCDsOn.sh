#!/bin/bash
#
#
#-------------------------------------
#
#- Python that works with mutils (may need updates)
#PATH=/gpfs/slac/lsst/fs2/u1/dh/software/centos7-gcc48/bin:${PATH}
#PATH=/gpfs/slac/lsst/fs2/u1/devel/marshall/cfitsio/bin:${PATH}
#PATH=/gpfs/slac/lsst/fs2/u1/devel/marshall/anaconda3-2023.03/bin:${PATH}
#PATH=/gpfs/slac/lsst/fs2/u1/devel/marshall/mutils/bin:${PATH}
#export PATH

function usage {
  cat <<-EOM
  Description: Plot CCD Power On's as they occur by monitoring the log file (arg1)
  Usage: ${0##*/} logfilepath
    logfilepath is full path to active logfile
    note must run on logfile host
  Options:
    -h (print help msg)
    -s (save the plots in /tmp as png)
EOM
exit 1
}
#-- process commandline options
#
duration=
while getopts "hs" Option
do
  case $Option in
    h  ) usage;;
    s  ) save="yes";;
    *  ) usage;;   # Default.
  esac
done
shift $((OPTIND - 1))

if [ $# -ne 1 ] ; then
    usage
fi

cmd="plotPowerCCDsOn.sh -w"
if [ $save ] ; then
    cmd="plotPowerCCDsOn.sh -s -w"
fi

tail -5f $1 |\
stdbuf -oL gawk '/Command\(TokenizedCommand\) powerCCDsOn/ {printf("%s %s -- powerCCDsOn\n",$10, $1);}' |\
stdbuf -oL sed 's/\(^.*\): \[\(.*\)\].*\(powerCCDsOn\).*/\1 \2/'
#xargs -L1 ${cmd}
#xargs -L1 ${cmd}  2>/dev/null
