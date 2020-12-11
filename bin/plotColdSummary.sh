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
    -t (use time averaged data, good on long queries)
    -s saveFileName
    -d <duration>  [default is 10m]
EOM
exit 1
}
#-- process commandline options
#
duration=
savefile=
while getopts "htd:s:" Option
do
  case $Option in
    h  ) usage;;
    d  ) duration=$OPTARG;;
    t  ) timebins="yes";;
    s  ) savefile=$OPTARG;;
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
regexes+=('^[rh].*[xg]/Cold1/(.*Cooling$|.*Heat$)')
regexes+=('^[rh].*[xg]/Cold2/(.*Cooling$|.*Heat$)')
#
regexes+=('^refrig/Cold1/(.*Tmp$|.*Tmp_M)')
regexes+=('^refrig/Cold2/(.*Tmp$|.*Tmp_M)')
#
regexes+=('^hex/Cold1/.*Tmp$')
regexes+=('^hex/Cold2/.*Tmp$')
#
regexes+=('^hex/Cold1/.*Prs')
regexes+=('^hex/Cold2/.*Prs')
#

if [ $duration"XXX" == "XXX" ] ; then
      duration=10m
fi

sv=
if [ $savefile"XXX" != "XXX" ] ; then
      sv="--save ${savefile}"
fi

if [ $timebins ] ; then
      timebins='--timebins'
fi

echo trender.py ${st} ${sv} --dur ${duration} ${timebins} --title "Cold Refrig Summary" --plot --layout 4x2 --outside --overlayregex -- "${regexes[@]}"

