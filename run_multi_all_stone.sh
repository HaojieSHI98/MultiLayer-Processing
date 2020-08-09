cores=$1

objects=$2

updateratio=$3

quratio=$4

updateratio2=$5

quratio2=$6

layer=$7

display=$8

configstrs=$9

echo $quratio

echo $cores

echo $objects


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


/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 10 -testtime 20 -testtime2 20 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${updateratio2}_${quratio2}_time_10_20 -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -query2 $querynum2 -insert2 $insertnum2 -delete2 $deletenum2 -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6 -configstr $9
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
