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
    -v (verbose output)
    -e (echo command only)
EOM
exit 1
}
#-- process commandline options
#
duration=
savefile=
while getopts "evhtd:s:" Option
do
  case $Option in
    h  ) usage;;
    d  ) duration=$OPTARG;;
    t  ) timebins="yes";;
    s  ) savefile=$OPTARG;;
    v  ) verbose="yes";;
    e  ) echocmd="yes";;
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

if [ $verbose ] ; then
    verbose='--debug'
fi


declare -a regexes
regexes+=('thermal/Cold.*/CLP-RTD-.5')
regexes+=('focal-plane/R.4/Reb[2G]/Temp2')
regexes+=('thermal/Cold.*/CLP-RTD-.[23]')
regexes+=('focal-plane/R[04][13]/Reb2/Temp2')
regexes+=('thermal/Cold.*/CLP-RTD-.0')
regexes+=('focal-plane/R.0/Reb[2G]/Temp2')
regexes+=('hex/Cold1/(EvapExitTmp|PreExpnTmp|EvapSuperHeat)')
regexes+=('hex/Cold2/(EvapExitTmp|PreExpnTmp|EvapSuperHeat)')
regexes+=('^[hr].*/Cold./(Return|Suction)Prs$')
regexes+=('^[hr].*/Cold./(Dis.*|Sup.*)Prs$')
regexes+=('thermal/Trim_Htrs/ColdHtr[035]_P')
regexes+=('refrig/Cold./EEPRValvePosn')

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

if [ $echocmd ] ; then
      echo trender.py ${verbose} ${st} ${sv} --dur ${duration} ${timebins} --title "Cold Summary" --plot --outside --overlayregex -- "${regexes[@]}"
      exit
fi

trender.py --layout 6x2 ${verbose} ${st} ${sv} --dur ${duration} ${timebins} --title "Cold Thermal Controls Summary" --plot --outside --overlayregex -- "${regexes[@]}"

