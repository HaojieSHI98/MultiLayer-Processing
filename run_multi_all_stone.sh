cores=$1

objects=$2

updateratio=$3

quratio=$4

layer=$5

query_thread=$6

update_thread=$7

echo $quratio

echo $cores

echo $objects

echo $query_thread

echo $update_thread

echo $query_cost

echo $insert_cost

echo $delete_cost

fullusecore=`expr "$cores" - "2"`

updatenum=`expr "${objects}" \* "${updateratio}"`

insertnum=`expr "$updatenum" / "2"`

deletenum=`expr "$updatenum" / "2"`

querynum=`expr "$updatenum" \* "$quratio"`
echo $cores
echo $querynum
echo $insertnum
echo $deletenum

/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 50 -testtime 10 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${query_thread}_${update_thread}_${query_cost}_${insert_cost}_${delete_cost} -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -init $objects -layer 1 -toaintype q -network NY -querythread $6 -updatethread $7
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -configtime 50 -testtime 10 -method dijk -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${updateratio}_${quratio}_${query_thread}_${update_thread}_${query_cost}_${insert_cost}_${delete_cost} -single_aggr 1 -query $querynum -insert $insertnum -delete $deletenum -init $objects -layer 1 -toaintype q -network NY -querythread $6 -updatethread $7
