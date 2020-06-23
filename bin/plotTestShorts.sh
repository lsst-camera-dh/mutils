#!/bin/bash
#
#------------------------------------------------------------------------

function usage {
  cat <<-EOM
  Usage ${0##*/} rebPath startTime
    rebPath ~ <subsystem>/<bay>/Reb[012]
    startTime ~ 2020-06-19T11:00:41-07:00
    quote time if it contains spaces
EOM
}

if [ $# -ne 2 ]; then
    usage
    exit 1
fi

rebPath=$1
st=$2

declare -a regexes
regexes[0]=${1}'/[CO].*I'       # board currents
regexes[1]=${1}'/[PSR].*[UL]$'  # clock levels
regexes[2]=${1}'/S../.*V$'      # bias voltages

trender.py --lay 3x1 --out --start "${st}" --title testCCDShorts:${1} --overlayreg --plot --dur 9s --fmt 'o-' -- ${regexes[@]}

# All channels for Reb0
#
# R22/Reb0/Temp1, R22/Reb0/Temp2, R22/Reb0/Temp3, R22/Reb0/Temp4, R22/Reb0/Temp5,
# R22/Reb0/Temp6, R22/Reb0/Temp7, R22/Reb0/Temp8, R22/Reb0/Temp9, R22/Reb0/Temp10,
# R22/Reb0/AspicU/Temp0, R22/Reb0/AspicL/Temp0,
# R22/Reb0/AspicU/Temp1, R22/Reb0/AspicL/Temp1,
# R22/Reb0/AspicU/Temp2, R22/Reb0/AspicL/Temp2,
# R22/Reb0/HVBiasSwitch,
# R22/Reb0/S00/Temp, R22/Reb0/S01/Temp, R22/Reb0/S02/Temp, R22/Reb0/RTDTemp,
# R22/Reb0/HtrV, R22/Reb0/HtrW,
# R22/Reb0/DigV, R22/Reb0/DigI, R22/Reb0/AnaV, R22/Reb0/AnaI,
# R22/Reb0/ClkHV, R22/Reb0/ClkHI, R22/Reb0/ClkLV, R22/Reb0/ClkLI,
# R22/Reb0/ODV, R22/Reb0/ODI, R22/Reb0/Power,
# R22/Reb0/PClkU, R22/Reb0/PClkL, R22/Reb0/SClkU, R22/Reb0/SClkL, R22/Reb0/RGU, R22/Reb0/RGL,
# R22/Reb0/RefP12, R22/Reb0/RefN12,
# R22/Reb0/S00/ODV, R22/Reb0/S00/OGV, R22/Reb0/S00/RDV, R22/Reb0/S00/GDV,
# R22/Reb0/S01/ODV, R22/Reb0/S01/OGV, R22/Reb0/S01/RDV, R22/Reb0/S01/GDV,
# R22/Reb0/S02/ODV, R22/Reb0/S02/OGV, R22/Reb0/S02/RDV, R22/Reb0/S02/GDV,
# R22/Reb0/Ref05V, R22/Reb0/Ref15V, R22/Reb0/Ref25V, R22/Reb0/Ref125V,
# R22/Reb0/S00/Seg00/I, R22/Reb0/S00/Seg01/I, R22/Reb0/S00/Seg02/I, R22/Reb0/S00/Seg03/I,
# R22/Reb0/S00/Seg04/I, R22/Reb0/S00/Seg05/I, R22/Reb0/S00/Seg06/I, R22/Reb0/S00/Seg07/I,
# R22/Reb0/S00/Seg10/I, R22/Reb0/S00/Seg11/I, R22/Reb0/S00/Seg12/I, R22/Reb0/S00/Seg13/I,
# R22/Reb0/S00/Seg14/I, R22/Reb0/S00/Seg15/I, R22/Reb0/S00/Seg16/I, R22/Reb0/S00/Seg17/I,
# R22/Reb0/S01/Seg00/I, R22/Reb0/S01/Seg01/I, R22/Reb0/S01/Seg02/I, R22/Reb0/S01/Seg03/I,
# R22/Reb0/S01/Seg04/I, R22/Reb0/S01/Seg05/I, R22/Reb0/S01/Seg06/I, R22/Reb0/S01/Seg07/I,
# R22/Reb0/S01/Seg10/I, R22/Reb0/S01/Seg11/I, R22/Reb0/S01/Seg12/I, R22/Reb0/S01/Seg13/I,
# R22/Reb0/S01/Seg14/I, R22/Reb0/S01/Seg15/I, R22/Reb0/S01/Seg16/I, R22/Reb0/S01/Seg17/I,
# R22/Reb0/S02/Seg00/I, R22/Reb0/S02/Seg01/I, R22/Reb0/S02/Seg02/I, R22/Reb0/S02/Seg03/I,
# R22/Reb0/S02/Seg04/I, R22/Reb0/S02/Seg05/I, R22/Reb0/S02/Seg06/I, R22/Reb0/S02/Seg07/I,
# R22/Reb0/S02/Seg10/I, R22/Reb0/S02/Seg11/I, R22/Reb0/S02/Seg12/I, R22/Reb0/S02/Seg13/I,
# R22/Reb0/S02/Seg14/I, R22/Reb0/S02/Seg15/I, R22/Reb0/S02/Seg16/I, R22/Reb0/S02/Seg17/I,
