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

if [ $duration"XXX" == "XXX" ] ; then
      duration=10m
fi

TRENDER=trender.py
OPTIONS="--site summit --dur ${duration} --stat --nosort -- "
GAWK=gawk
GOPTIONS='{printf("%7.3g %-25.25s %3d/%-3d %5.5s\n", $3, $8, ($3 ? 100.0*$6/$3 : 0.0), ($3 ? 100.0*$5/$3 : 0.0), $9);}'

printf "#%6s %-25.25s %7.7s %5.5s\n" "median" "path" "+/-%" "units"
printf "#\n"

# vacuum
printf "# vacuum\n"
declare -a regexes
regexes+=('vacuum/(Cryo|HX)/[CH].*Vac')
regexes+=('vacuum/Cryo/[TF].*Vac')
regexes+=('vacuum/HX/[TF].*Vac')
regexes+=('vacuum/Cryo/TurboP.*')
regexes+=('vacuum/HX/TurboP.*')
regexes+=('vacuum/Inst/PumpCartPressure')
$TRENDER ${OPTIONS} "${regexes[@]}" | $GAWK "/^#/ {next;}; ${GOPTIONS}"
unset regexes
printf "#\n"

# refrig 
printf "# refrig\n"
declare -a regexes
regexes+=('chiller1/Chiller/FluidTemperature')
regexes+=('chiller1/Chiller/CoolPercentage')
regexes+=('chiller1/Maq20/Stg2DeSuHtrOut')
regexes+=('chiller1/Maq20/Glyc(Chiller.*|InputFlow)')
regexes+=('^refrig/Cryo.*/DischrgPrs')
regexes+=('^hex/Cryo.*/HexRtrnTmp')
$TRENDER ${OPTIONS} "${regexes[@]}" | $GAWK "/^#/ {next;}; ${GOPTIONS}"
unset regexes
printf "#\n"

# Body
printf "# Body\n"
declare -a regexes
regexes+=('utiltrunk/Body/Shtr.*AirVel')
regexes+=('utiltrunk/Body/(A|Dome|Shtr|VPP|Chgr|L2XM|CamHousXM).*emp')
$TRENDER ${OPTIONS} "${regexes[@]}" | $GAWK "/^#/ {next;}; ${GOPTIONS}"
unset regexes
printf "#\n"

# thermal
printf "# thermal\n"
declare -a regexes
regexes+=('focal-plane/R[024][13]/Reb./S[1W][1]*/Temp')
regexes+=('thermal/Cold_Temp/AvgColdTemp')
regexes+=('thermal/Cryo_Temp/AvgCryoTemp')
regexes+=('thermal/Trim_Htrs/.*Total_P')
regexes+=('focal-plane/RebTotalPower')
$TRENDER ${OPTIONS} "${regexes[@]}" | $GAWK "/^#/ {next;}; ${GOPTIONS}"
unset regexes
printf "#\n"

# UT
printf "# UT\n"
declare -a regexes
regexes+=('quadbox/PDU_[24].*/.*_T$')
regexes+=('^utiltrunk/UT/([AF].*Temp|Cool[FP].*|FanSpeed)')
$TRENDER ${OPTIONS} "${regexes[@]}" | $GAWK "/^#/ {next;}; ${GOPTIONS}"
unset regexes
printf "#\n"


