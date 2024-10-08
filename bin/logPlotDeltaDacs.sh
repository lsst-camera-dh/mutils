#!/bin/bash
#
#
#-------------------------------------
#
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

cmd="logPlotLoadDeltaDacs.sh -w"
if [ $save ] ; then
    cmd="logPlotLoadDeltaDacs.sh -s -w"
fi

tail -5f $1 |\
stdbuf -oL gawk '/Command\(TokenizedCommand\) loadDeltaDacs/ {printf("%s %s -- loadDeltaDacs\n",$9, $1);}' |\
stdbuf -oL sed 's/\(^.*\): \[\(.*\)\].*\(loadDeltaDacs\).*/\1 \2/' |\
xargs -L1 ${cmd}  2>/dev/null
