#!/bin/bash
#
#------------------------------------------------------------------------
function usage {
  cat <<-EOM
  Usage ${0##*/} rebPath startTime
    rebPath ~ <subsystem>/<bay>/Reb[012]
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
while getopts "hswd:" Option
do
  case $Option in
    h  ) usage;;
    d  ) duration=$OPTARG;;
    s  ) savePlot="yes";;
    w  ) waitTime="yes";;
    *  ) ${ECHO} "Unimplemented option chosen.";;   # Default.
  esac
done
shift $((OPTIND - 1))

declare -a regexes
regexes+=(${1}'/[DA].*I')      # digital & analog currents
regexes+=(${1}'/[CO].*I')       # clocks and OD currents
regexes+=(${1}'/[PSR].*[UL]$')  # clock levels
regexes+=(${1}'/S.*/.*V$')      # bias voltages

if [ $duration"XXX" == "XXX" ] ; then
      duration=8s
fi

if [ $waitTime ] ; then
      sleep $duration
      sleep 8
fi

sarg=
if [ $savePlot ] ; then
    sarg=" --save "$(echo -n $1 | sed 's?/?_?g' | sed 's?\^??' )"_TestShorts.png"
fi

site=
if [[ $HOSTNAME"XXX" =~ .*slac.stanford.edu"XXX" ]] ; then
    site="slac"
else
    site="localhost"
fi

trender.py --site ${site} --lay 4x1 --out ${sarg} --start "${2}" --title "testCCDShorts:${1}" --overlayreg --plot --dur $duration --fmt 'o-' -- "${regexes[@]}"
