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
    -e (echo command)
    -d <duration>  [default is 10m]
EOM
exit 1
}
#-- process commandline options
#
duration=
savefile=
while getopts "hted:s:" Option
do
  case $Option in
    h  ) usage;;
    d  ) duration=$OPTARG;;
    t  ) timebins="yes";;
    s  ) savefile=$OPTARG;;
    e  ) doecho="yes";;
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
regexes+=('^hex/Cold1/.*Tmp')
regexes+=('^hex/Cold2/.*Tmp')
regexes+=('^thermal/Cold_Temp/CLP-RTD-..')
regexes+=('^thermal/Cryo_Temp/CYP-RTD-.[24]')
regexes+=('^refrig/Cold1/([^D].*Tmp_M$|[^C].*Tmp$)')
regexes+=('^refrig/Cold2/([^D].*Tmp_M$|[^C].*Tmp$)')
regexes+=('^[hr].*[xg]/Cold./(Supply|Dischrg).*Prs')
regexes+=('^[hr].*[xg]/Cold./(Suction|Return).*Prs')
regexes+=('^refrig/Cold./CompPower$')
regexes+=('^[tf].*[le]/.*(ColdTotal_P$|CryoTotal_P$|RebTotalPower$)')
regexes+=('^refrig/Cold./.*ValvePosn')
regexes+=('^refrig/Cold./CoolFlowRate')

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

if [ $doecho ] ; then
    echo trender.py ${st} ${sv} --dur ${duration} ${timebins} --title "Refrig Summary" --plot --layout 8x2 --outside --overlayregex -- "${regexes[@]}"
    exit
fi

trender.py ${st} ${sv} --dur ${duration} ${timebins} --title "Refrig Summary" --plot --layout 8x2 --outside --overlayregex -- "${regexes[@]}"

