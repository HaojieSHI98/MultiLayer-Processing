cores=$1

objects=$2

updateratio=$3

quratio=$4

layer=$5

display=$6

echo $quratio

echo $cores

echo $objects


fullusecore=`expr "$cores" - "2"`

updatenum=`expr "${objects}" \* "${updateratio}"`

insertnum=`expr "$updatenum" / "2"`

deletenum=`expr "$updatenum" / "2"`

querynum=`expr "$updatenum" \* "$quratio"`
echo $cores
echo $querynum
echo $insertnum
echo $deletenum

/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 50 -testtime 10 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio} -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 50 -testtime 10 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${query_thread}_${update_thread}_${query_cost}_${insert_cost}_${delete_cost} -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -init $objects -layer 1 -toaintype q -network NY -querythread $6 -updatethread $7
