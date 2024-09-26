#!/bin/bash
#
#------------------------------------------------------------------------
function usage {
  cat <<-EOM
  Usage ${0##*/} <raftpath> [stopTime]
    stopTime ~ 2020-06-19T11:00:41-07:00
    quote time if it contains spaces
  Options:
    -h (print help msg)
    -d <duration>  [default is 10m]
EOM
exit 1
}
#-- process commandline options
#
duration=
while getopts "hd:" Option
do
  case $Option in
    h  ) usage;;
    d  ) duration=$OPTARG;;
    *  ) ${ECHO} "Unimplemented option chosen.";;   # Default.
  esac
done
shift $((OPTIND - 1))

if [ $# -gt 2 ]; then
    usage
fi

if [ $2"XXX" == "XXX" ] ; then
    st="--stop "$(date --iso-8601=s)
else
    st="--stop ${2}"
fi

declare -a regexes
regexes[0]=${1}'/Reb./[DACO].*V$'  # board voltages
regexes[1]=${1}'/Reb./[DACO].*I$'  # board currents
regexes[2]=${1}'/Reb./[PSR].*[UL]$'  # clock levels
regexes[3]=${1}'/Reb./S.*/.*V$'      # bias voltages
regexes[4]=${1}'/Reb./Temp[0-9]*$' # reb temps
regexes[5]=${1}'/Reb./Aspic[UL]/Temp.$' # aspic temps
regexes[6]=${1}'/Reb./S../Temp$'
regexes[7]=${1}'/Reb./(HtrW|Power)$'

if [ $duration"XXX" == "XXX" ] ; then
      duration=10m
fi

trender.py --deb ${st} --dpi 200 --dur $duration --title "Raft Summary: ${1}" --plot --layout 4x2 --outside --overlayregex -- "${regexes[@]}"
