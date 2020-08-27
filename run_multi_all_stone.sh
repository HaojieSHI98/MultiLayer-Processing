cores=$1

objects=$2

configstrs=$3

teststrs=$4

layer=$5

display=$6

mode=$7

echo $quratio

echo $cores

echo $objects

echo $configstr


fullusecore=`expr "$cores" - "2"`

updatenum=`expr "${objects}" \* "${updateratio}"`

insertnum=`expr "$updatenum" / "2"`

deletenum=`expr "$updatenum" / "2"`

querynum=`expr "$updatenum" \* "$quratio"`

updatenum2=`expr "${objects}" \* "${updateratio2}"`

insertnum2=`expr "$updatenum2" / "2"`

deletenum2=`expr "$updatenum2" / "2"`

querynum2=`expr "$updatenum2" \* "$quratio2"`

echo $cores
echo $querynum
echo $insertnum
echo $deletenum
echo $querynum2
echo $insertnum2
echo $deletenum2


/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -teststr $4  -method toain -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${configstr} -single_aggr 1  -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6 -configstr $3 -mode $7
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6

#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 50 -testtime 10 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${query_thread}_${update_thread}_${query_cost}_${insert_cost}_${delete_cost} -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -init $objects -layer 1 -toaintype q -network NY -querythread $6 -updatethread $7
