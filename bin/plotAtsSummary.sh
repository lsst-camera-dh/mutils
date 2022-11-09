#!/bin/bash
#
#------------------------------------------------------------------------
function usage {
  cat <<-EOM
  Usage ${0##*/} <rebpath> [stopTime]
    stopTime ~ 2020-06-19T11:00:41-07:00
    quote time if it contains spaces
  Options:
    -h (print help msg)
    -d <duration>  [default is 10m]
    -t (use time averaged data, good on long queries)
EOM
exit 1
}
#-- process commandline options
#
duration=
while getopts "hd:t" Option
do
  case $Option in
    h  ) usage;;
    d  ) duration=$OPTARG;;
    t  ) timebins="yes";;
    *  ) ${ECHO} "Unimplemented option chosen.";;   # Default.
  esac
done
shift $((OPTIND - 1))

if [ $# -gt 2 ] ; then
    usage
elif [ $# -eq 0 ] ; then
    usage
fi

if [ $2"XXX" == "XXX" ] ; then
    st="--stop "$(date --iso-8601=s)
else
    st="--stop ${2}"
fi

declare -a regexes
regexes[0]=${1}'/[DACO].*V$'  # board voltages
regexes[1]=${1}'/[DACO].*I$'  # board currents
regexes[2]=${1}'/[PSR].*[UL]$'  # clock levels
regexes[3]=${1}'/S.*/.*V$'      # bias voltages
regexes[4]=${1}'/Temp[0-9]*$' # reb temps
#regexes[5]=${1}'/Aspic[UL]/Temp.$' # aspic temps

if [ $duration"XXX" == "XXX" ] ; then
      duration=10m
fi

if [ $timebins ] ; then
      timebins='--timebins'
fi

#echo trender.py --site ats ${st} --dur $duration ${timebins} --title \"Reb Summary: ${1}\" --plot --layout 3x2 --outside --overlayregex -- \"${regexes[@]}\"
trender.py --site ats ${st} --dur $duration ${timebins} --title "Reb Summary: ${1}" --plot --layout 3x2 --outside --overlayregex -- "${regexes[@]}"

