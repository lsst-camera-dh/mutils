#!/bin/bash
#
#------------------------------------------------------------------------
function usage {
  cat <<-EOM
  Usage ${0##*/} [startTime]
    startTime ~ 2020-06-19T11:00:41-07:00
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

if [ $# -gt 1 ]; then
    usage
fi

if [ $1"XXX" == "XXX" ] ; then
    st="--stop "$(date --iso-8601=s)
else
    st="--start ${1}"
fi

declare -a regexes
regexes[0]='^[rh].*[xg]/Cold1/(.*Cooling$|.*Heat$)' 
regexes[1]='^[rh].*[xg]/Cold2/(.*Cooling$|.*Heat$)'
regexes[2]='^refrig/Cold1/(.*Tmp_M$|[^A^C].*Tmp$)'
regexes[3]='^refrig/Cold2/(.*Tmp_M$|[^A^C].*Tmp$)' 
regexes[4]='^hex/Cold1/.*Tmp'
regexes[5]='^hex/Cold2/.*Tmp'
regexes[6]='^[hr].*[xg]/Cold./(Supply|Dischrg).*Prs'
regexes[7]='^[hr].*[xg]/Cold./(Suction|Return).*Prs'
regexes[8]='^[rf].*/Cold./CompPower$'
regexes[9]='^[tf].*[es]/(ColdTotal_P$|CryoTotal_P$|RebTotalPower$)'
regexes[10]='^thermal/Cold_Temp/CLP-RTD-..'
regexes[11]='^thermal/Cryo_Temp/CYP-RTD-.[24]'


if [ $duration"XXX" == "XXX" ] ; then
      duration=10m
fi

trender.py ${st} --dur $duration --title "Refrig Summary" --plot --layout 6x2 --outside --overlayregex -- "${regexes[@]}"
