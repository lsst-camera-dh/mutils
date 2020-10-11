#!/bin/bash
#
#------------------------------------------------------------------------
function usage {
  cat <<-EOM
  Usage ${0##*/} rebPath startTime
    rebPath ~ <subsystem>/<bay>/Reb[012]
    startTime ~ 2020-06-19T11:00:41-07:00
    quote time if it contains spaces
  Options:
    -h (print help msg)
    -d <duration>
EOM
exit 1
}
if [ $# -lt 2 ]; then
    usage
fi

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

declare -a regexes
regexes[0]=${1}'/[CO].*I'       # board currents
regexes[1]=${1}'/[PSR].*[UL]$'  # clock levels
regexes[2]=${1}'/S../.*V$'      # bias voltages
if [ $duration"XXX" == "XXX" ] ; then
      duration=9s
fi

trender.py --lay 3x1 --out --start "${2}" --title "testCCDShorts:${1}" --overlayreg --plot --dur $duration --fmt 'o-' -- "${regexes[@]}"
