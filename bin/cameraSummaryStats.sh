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
savefile=
while getopts "htd:s:c:p:" Option
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
    st="--stop ${1}"
fi

declare -a regexes
regexes+=('thermal/Cold_Temp/AvgColdTemp')
regexes+=('thermal/Cryo_Temp/AvgCryoTemp')
regexes+=('focal-plane/R[024][13]/Reb./S[1W][1]*/Temp')
regexes+=('vacuum/Cryo/[CTF].*Vac')
regexes+=('vacuum/Cryo/TurboP.*')
regexes+=('vacuum/HX/[HTF].*Vac')
regexes+=('vacuum/HX/TurboP.*')
regexes+=('vacuum/Inst/PumpCartPressure')
regexes+=('quadbox/PDU_[24].*/.*_T$')
regexes+=('chiller/Chiller/FluidTemperature')
regexes+=('chiller/Chiller/CoolPercentage')
regexes+=('chiller/Maq20/Stg2DeSuHtrOut')
regexes+=('focal-plane/RebTotalPower')
regexes+=('thermal/Trim_Htrs/.*Total_P')
regexes+=('^utiltrunk/UT/([AF].*Temp|Cool.*|FanSpeed)')
regexes+=('^refrig/Cryo.*/DischrgPrs')
regexes+=('^hex/Cryo.*/HexRtrnTmp')

if [ $duration"XXX" == "XXX" ] ; then
      duration=10m
fi

gawk 'BEGIN {printf("#%6s %-25.25s %7.7s %5.5s\n", "median", "path", "+/-%", "units");}'
trender.py --site summit --dur ${duration} --stat -- "${regexes[@]}" |\
    gawk '/^#/ {next;}; {printf("%7.3g %-25.25s %3d/%-3d %5.5s\n", $3, $8, ($3 ? 100.0*$6/$3 : 0.0), ($3 ? 100.0*$5/$3 : 0.0), $9);}' |\
    sort -r -k 2


