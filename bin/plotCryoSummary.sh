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
    -p <dpi>
    -c cryos (regex eg [02], default is . for all)
EOM
exit 1
}
#-- process commandline options
#
duration=
savefile=
while getopts "htd:s:c:p:" Option
do
  case $Option in
    h  ) usage;;
    d  ) duration=$OPTARG;;
    c  ) cryos=$OPTARG;;
    t  ) timebins="yes";;
    s  ) savefile=$OPTARG;;
    p  ) dpi=$OPTARG;;
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

if [ $cryos"XXX" == "XXX" ] ; then
      cryos='.'
fi

declare -a regexes
regexes+=('^refrig/Cryo'${cryos}'/CompPower')
regexes+=('^refrig/Cryo'${cryos}'/SuctionPrs')
regexes+=('^refrig/Cryo'${cryos}'/DischrgPrs')
regexes+=('^hex/Cryo'${cryos}'/.*C3.*Tmp')
regexes+=('^hex/Cryo'${cryos}'/.*C4.*Tmp')
regexes+=('^hex/Cryo'${cryos}'/EvapExitTmp')
regexes+=('^hex/Cryo'${cryos}'/HexRtrnTmp')
regexes+=('^(thermal/Cryo_Temp/CYP-RTD-(02|12|42|52)|fo.*/R[0134]2/Reb1/S11/Temp)')
#regexes+=('^[tf].*[le]/.*(C.*Total_P$|RebTotalPower$)')

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

dpis=
if [ $dpi"XXX" != "XXX" ] ; then
      dpis="--dpi ${dpi}"
fi

trender.py ${st} ${sv} ${dpis} --dur ${duration} ${timebins} --title "CryoRefrig Summary" --plot --outside --overlayregex -- "${regexes[@]}"

