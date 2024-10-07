#!/bin/csh

set qfile = $1
set rfile = ${qfile}.raw
set lfile = ${qfile}.l0
set station = $2

ql02raw_seawifs <<EOF
$qfile
$rfile
EOF

raw2l0_seawifs < $rfile > $lfile

set pref = `l0info_seawifs $lfile | grep Start_Zulu | cut -b 26-38`
set ofile = "S${pref}.L0_${station}"

mv $lfile $ofile
rm $rfile

exit 0
