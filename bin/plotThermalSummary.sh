#!/bin/bash
#
#------------------------------------------------------------------------
function usage {
  cat <<-EOM
  Plot CCD temps for all rafts
  Usage ${0##*/} [stopTime]
    stopTime ~ 2020-06-19T11:00:41-07:00
    quote time if it contains spaces
  Options:
    -h (print help msg)
    -d <duration>  [default is 10m]
    -t (use time averaged data, good on long queries)
    -p (only plot temps used in PID controls)
    -e echo command line only
EOM
exit 1
}
#-- process commandline options
#
duration=
while getopts "ehd:tp" Option
do
  case $Option in
    h  ) usage;;
    e  ) echocmd="yes";;
    d  ) duration=$OPTARG;;
    t  ) timebins="yes";;
    p  ) pidonly="yes";;
    *  ) ${ECHO} "Unimplemented option chosen.";;   # Default.
  esac
done
shift $((OPTIND - 1))

if [ $# -gt 1 ] ; then
    usage
fi

if [ $1"XXX" == "XXX" ] ; then
    st="--stop "$(date --iso-8601=s)
else
    st="--stop ${1}"
fi

declare -a regexes
regexes+=('focal-plane/R../Reb[1W]/S.1/Temp')
regexes+=('focal-plane/R../Reb[02W]/HtrW')
regexes+=('focal-plane/R../Reb1/Temp2$')
regexes+=('focal-plane/R../Reb1/(Temp6$|Aspic./Temp.)')
regexes+=('^thermal/(Cryo_Temp/CYP-RTD-.[24]$|Grid_Temp/GFX-RTD-0[13]$|Grid_Temp/GRD-RTD-0.$)')
regexes+=('^thermal/Cold_Temp/CLP-RTD-..')
regexes+=('^thermal/Grid_Temp/GFX-RTD-0[24]')
regexes+=('^thermal/Trim_Htrs/C.*Total_P')

#badrtds:
declare -a badCCDrtds
badCCDrtds+=('focal-plane/R11/Reb0/S02/Temp')
badCCDrtds+=('focal-plane/R41/Reb2/S20/Temp')
badCCDrtds+=('focal-plane/R41/Reb2/S21/Temp')
badCCDrtds+=('focal-plane/R41/Reb2/S22/Temp')
#nonPIDrtds
declare -a nonPIDrtds
nonPIDrtds+=('focal-plane/R02/Reb0/S02/Temp')
nonPIDrtds+=('focal-plane/R10/Reb1/S12/Temp')
nonPIDrtds+=('focal-plane/R12/Reb2/S20/Temp')
nonPIDrtds+=('focal-plane/R12/Reb2/S22/Temp')
nonPIDrtds+=('focal-plane/R13/Reb2/S21/Temp')
nonPIDrtds+=('focal-plane/R21/Reb0/S00/Temp')
nonPIDrtds+=('focal-plane/R21/Reb2/S21/Temp')
nonPIDrtds+=('focal-plane/R24/Reb1/S11/Temp')
nonPIDrtds+=('focal-plane/R24/Reb2/S22/Temp')
nonPIDrtds+=('focal-plane/R31/Reb2/S21/Temp')
nonPIDrtds+=('focal-plane/R31/Reb2/S22/Temp')
nonPIDrtds+=('focal-plane/R33/Reb1/S11/Temp')
nonPIDrtds+=('focal-plane/R33/Reb2/S22/Temp')
nonPIDrtds+=('focal-plane/R34/Reb0/S01/Temp')
nonPIDrtds+=('focal-plane/R34/Reb0/S02/Temp')
nonPIDrtds+=('focal-plane/R41/Reb0/S01/Temp')
nonPIDrtds+=('focal-plane/R42/Reb0/S02/Temp')
nonPIDrtds+=('focal-plane/R43/Reb0/S01/Temp')
nonPIDrtds+=('focal-plane/R43/Reb1/S12/Temp')
nonPIDrtds+=('focal-plane/R43/Reb2/S21/Temp')

if [ $duration"XXX" == "XXX" ] ; then
      duration=10m
fi

if [ $timebins ] ; then
      timebins='--timebins'
fi

if [ $pidonly ] ; then
    badCCDrtds+=( ${nonPIDrtds[@]} )
fi

if [ $echocmd ] ; then
    echo trender.py ${st} --dur $duration ${timebins} --title \"CryoStat Thermal Summary\" --layout 4x2 --plot --outside --reject \"${badCCDrtds[@]}\" --overlayregex -- \"${regexes[@]}\"
    exit 0
fi

trender.py ${st} --dur $duration ${timebins} --title "CryoStat Thermal Summary" --layout 4x2 --plot --outside --reject "${badCCDrtds[@]}" --overlayregex -- "${regexes[@]}"
