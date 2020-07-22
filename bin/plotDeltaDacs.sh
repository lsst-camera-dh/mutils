#!/bin/bash
#
#------------------------------------------------------------------------

function usage {
  cat <<-EOM
  Usage ${0##*/} rebPath startTime
    rebPath ~ <subsystem>/<bay>/Reb[012]
    startTime ~ 2020-06-19T11:00:41-07:00
    quote time if it contains spaces
EOM
}

if [ $# -ne 2 ]; then
    usage
    exit 1
fi

declare -a regexes
regexes[0]=${1}'/[CO].*I'       # board currents
regexes[1]=${1}'/[PSR].*[UL]$'  # clock levels
regexes[2]=${1}'/S.*/.*V$'      # bias voltages

trender.py --lay 3x1 --out --start "${2}" --title loadDeltaDacs:${1} --overlayreg --plot --dur 6s --fmt "o-" -- ${regexes[@]}
