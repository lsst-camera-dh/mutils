#!/bin/bash
#
#------------------------------------------------------------------------
function usage {
  cat <<-EOM
  Usage ${0##*/} [stopTime]
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
while getopts "htd:" Option
do
  case $Option in
    h  ) usage;;
    d  ) duration=$OPTARG;;
    t  ) timebins="yes";;
    *  ) ${ECHO} "Unimplemented option chosen.";;   # Default.
  esac
done
shift $((OPTIND - 1))

if [ $# -gt 1 ]; then
    usage
fi

if [ $1"XXX" == "XXX" ] ; then
    st="--stop "$(date --iso-8601=s)
else
    st="--stop ${1}"
fi

declare -a regexes
regexes+=('hex/Cryo1/.*Tmp')
regexes+=('hex/Cryo2/.*Tmp')
regexes+=('hex/Cryo3/.*Tmp')
regexes+=('hex/Cryo4/.*Tmp')
regexes+=('hex/Cryo5/.*Tmp')
regexes+=('hex/Cryo6/.*Tmp')

if [ $duration"XXX" == "XXX" ] ; then
      duration=10m
fi

if [ $timebins ] ; then
      timebins='--timebins'
fi

trender.py ${st} --dur ${duration} ${timebins} --title "Cryo Summary" --plot --layout 3x2 --outside --overlayregex -- "${regexes[@]}"
