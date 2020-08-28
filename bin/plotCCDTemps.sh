#!/bin/bash
#
#------------------------------------------------------------------------
function usage {
  cat <<-EOM
  Plot CCD temps for all rafts
  Usage ${0##*/} [stopTime]
    stopTime ~ 2020-06-19T11:00:41-07:00
    quote time if it contains spaces
  Options:
    -h (print help msg)
    -d <duration>  [default is 10m]
    -t (use time averaged data, good on long queries)
    -p (only plot temps used in PID controls)
EOM
exit 1
}
#-- process commandline options
#
duration=
while getopts "hd:tp" Option
do
  case $Option in
    h  ) usage;;
    d  ) duration=$OPTARG;;
    t  ) timebins="yes";;
    p  ) pidonly="yes";;
    *  ) ${ECHO} "Unimplemented option chosen.";;   # Default.
  esac
done
shift $((OPTIND - 1))

if [ $# -gt 1 ] ; then
    usage
fi

if [ $1"XXX" == "XXX" ] ; then
    st="--stop "$(date --iso-8601=s)
else
    st="--stop ${1}"
fi

declare -a regexes
regexes[0]='focal-plane/R../Reb./S../Temp'
#regexes[1]=${1}'/[DACO].*I$'  # board currents

#badrtds:
declare -a badCCDrtds
badCCDrtds[0]='focal-plane/R11/Reb0/S02/Temp'
badCCDrtds[1]='focal-plane/R41/Reb2/S20/Temp'
badCCDrtds[2]='focal-plane/R41/Reb2/S21/Temp'
badCCDrtds[3]='focal-plane/R41/Reb2/S22/Temp'

#nonPIDrtds
nonPIDrtds[0]='focal-plane/R12/Reb2/S20/Temp'
nonPIDrtds[1]='focal-plane/R12/Reb2/S22/Temp'
nonPIDrtds[2]='focal-plane/R13/Reb2/S21/Temp'
nonPIDrtds[3]='focal-plane/R21/Reb0/S00/Temp'
nonPIDrtds[4]='focal-plane/R21/Reb2/S21/Temp'
nonPIDrtds[5]='focal-plane/R24/Reb1/S11/Temp'
nonPIDrtds[6]='focal-plane/R24/Reb2/S22/Temp'
nonPIDrtds[7]='focal-plane/R31/Reb2/S21/Temp'
nonPIDrtds[8]='focal-plane/R31/Reb2/S22/Temp'
nonPIDrtds[9]='focal-plane/R33/Reb1/S11/Temp'
nonPIDrtds[10]='focal-plane/R33/Reb2/S22/Temp'
nonPIDrtds[11]='focal-plane/R34/Reb0/S01/Temp'
nonPIDrtds[12]='focal-plane/R34/Reb0/S02/Temp'
nonPIDrtds[13]='focal-plane/R41/Reb0/S01/Temp'
nonPIDrtds[14]='focal-plane/R43/Reb0/S01/Temp'
nonPIDrtds[15]='focal-plane/R43/Reb1/S12/Temp'
nonPIDrtds[16]='focal-plane/R43/Reb2/S21/Temp'
nonPIDrtds[17]='focal-plane/R02/Reb0/S02/Temp'
nonPIDrtds[18]='focal-plane/R42/Reb0/S02/Temp'
nonPIDrtds[19]='focal-plane/R10/Reb1/S12/Temp'

if [ $duration"XXX" == "XXX" ] ; then
      duration=10m
fi

if [ $timebins ] ; then
      timebins='--timebins'
fi

if [ $pidonly ] ; then
    badCCDrtds+=( ${nonPIDrtds[@]} )
fi

trender.py ${st} --dur $duration ${timebins} --title "Focal Plane CCD Temps" --plot --outside --reject "${badCCDrtds[@]}" --overlayregex -- "${regexes[@]}"


