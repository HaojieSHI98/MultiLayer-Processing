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


/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -teststr $4  -method toain -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${configstr} -single_aggr 1  -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6 -configstr $3 -mode $7 -updatethread 1
/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -teststr $4  -method toain -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${configstr} -single_aggr 1  -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6 -configstr $3 -mode $7 -updatethread 3
/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -teststr $4  -method toain -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${configstr} -single_aggr 1  -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6 -configstr $3 -mode $7 -updatethread 18
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -teststr $4  -method toain -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${configstr} -single_aggr 1  -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6 -configstr $3 -mode 0 -updatethread 1
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -teststr $4  -method toain -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${configstr} -single_aggr 1  -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6 -configstr $3 -mode 0 -updatethread 1
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -teststr $4  -method toain -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${configstr} -single_aggr 1  -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6 -configstr $3 -mode 1 -updatethread 1
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -teststr $4  -method toain -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${configstr} -single_aggr 1  -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6 -configstr $3 -mode 1 -updatethread 1
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -teststr $4  -method toain -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${configstr} -single_aggr 1  -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6 -configstr $3 -mode 1 -updatethread 1
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -teststr $4  -method toain -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${configstr} -single_aggr 1  -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6 -configstr $3 -mode 1 -updatethread 1
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -teststr $4  -method toain -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${configstr} -single_aggr 1  -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6 -configstr $3 -mode 1 -updatethread 1
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -teststr $4  -method toain -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${configstr} -single_aggr 1  -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6 -configstr $3 -mode $7 -updatethread 3
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -teststr $4  -method toain -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${configstr} -single_aggr 1  -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6 -configstr $3 -mode $7 -updatethread 3
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -teststr $4  -method toain -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${configstr} -single_aggr 1  -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6 -configstr $3 -mode $7 -updatethread 3
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -teststr $4  -method toain -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${configstr} -single_aggr 1  -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6 -configstr $3 -mode $7 -updatethread 3
#/home/siqiang/MPR/TOAIN/TOAIN -multicore $1 -part 1 -teststr $4  -method toain -threshold 0 -parmethod rand -out _NY_${cores}_${objects}_${configstr} -single_aggr 1  -init $objects -layer 1 -toaintype q -network NY -DISPLAY $6 -configstr $3 -mode $7 -updatethread 3