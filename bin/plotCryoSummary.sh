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

#channels:
#refrig/Cryo1/FanSpeed
#refrig/Cryo1/OilSepTmp
#refrig/Cryo1/SurgeTankTmp
#refrig/Cryo1/CompVoltage
#refrig/Cryo1/CompCurrent
#refrig/Cryo1/DischrgTmp_P
#refrig/Cryo1/SuctionTmp_P
#refrig/Cryo1/SuctionPrs
#refrig/Cryo1/OilLevel
#refrig/Cryo1/DischrgPrs
#refrig/Cryo1/CompPower
#refrig/Cryo1/DischrgTmp_M
#refrig/Cryo1/PhaseSepTmp
#refrig/Cryo1/SuctionTmp_M
#refrig/Cryo1/WaterInTmp
#refrig/Cryo1/WaterOutTmp
#channels:
#hex/Cryo1/C3ExitTmp
#hex/Cryo1/C4ExitTmp
#hex/Cryo1/LiquidPrs
#hex/Cryo1/PreC4Tmp
#hex/Cryo1/ReturnPrs
#hex/Cryo1/VaporPrs

cryos='[1356]'
#cryos='[12356]'
declare -a regexes
regexes+=('^refrig/Cryo'${cryos}'/CompPower')
regexes+=('^refrig/Cryo'${cryos}'/SuctionPrs')
regexes+=('^refrig/Cryo'${cryos}'/DischrgPrs')
regexes+=('^hex/Cryo'${cryos}'/.*C3.*Tmp')
regexes+=('^hex/Cryo'${cryos}'/EvapExitTmp')
regexes+=('^hex/Cryo'${cryos}'/HexRtrnTmp')
regexes+=('^(thermal/Cryo_Temp/CYP-RTD-(14|43|12|31)|fo.*/R2[0134]/Reb1/S11/Temp)')
regexes+=('^[tf].*[le]/.*(CryoTotal_P$|RebTotalPower$)')

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

trender.py ${st} ${sv} --dur ${duration} ${timebins} --title "CryoRefrig Summary" --plot --layout 4x2 --outside --overlayregex -- "${regexes[@]}"

