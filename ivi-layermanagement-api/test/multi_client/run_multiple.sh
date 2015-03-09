./ivi-layermanagement-api-test-multi --l300 --s200 --n100 --t100 $* > result.txt 2>&1 &
./ivi-layermanagement-api-test-multi --l700 --s500 --n100 --t100 $* > result2.txt 2>&1 &
./ivi-layermanagement-api-test-multi --l800 --s600 --n100 --t100 $* > result3.txt 2>&1
