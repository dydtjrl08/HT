#!/bin/bash

set -euo pipefail


if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <particle>"
  echo "Example: $0 proton"
  exit 1
fi


particle="$1"
radian=$2
echo "Starting all jobs in parallel for particle: $particle ..."

for dl in {5..50..5}; do
	 root -l -b -q "CurvatureHT_ErrorAnalRad.C(\"$particle\", $dl,$radian)" &
	 #root -l -b -q "CurvatureHT_ErrorAnalRad.C+(\"$particle\", $dl)" &
done

wait

echo "All jobs have finished."

#root -l -b -q 'CurvatureHT_ErrorAnalRad.C(5)' &
#root -l -b -q 'CurvatureHT_ErrorAnalRad.C(10)' &
#root -l -b -q 'CurvatureHT_ErrorAnalRad.C(15)' &
#root -l -b -q 'CurvatureHT_ErrorAnalRad.C(20)' &
#root -l -b -q 'CurvatureHT_ErrorAnalRad.C(25)' &
#root -l -b -q 'CurvatureHT_ErrorAnalRad.C(30)' &
#root -l -b -q 'CurvatureHT_ErrorAnalRad.C(35)' &
#root -l -b -q 'CurvatureHT_ErrorAnalRad.C(40)' &
#root -l -b -q 'CurvatureHT_ErrorAnalRad.C(45)' &
#root -l -b -q 'CurvatureHT_ErrorAnalRad.C(50)' &


#wait


#echo "All jobs have finished."



