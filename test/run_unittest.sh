#!/usr/bin/env bash

run_one_unittest() {
    if [ $# == 1 ];then
        fn=$1
        if [ -f $fn ];then
            echo "Running $fn"
            chmod +x $fn
            $fn
            status=$?
            if [ ! $status == 0 ];then
                exit
            fi
            echo ""
        else
            echo "!!WARN: $fn not found. Skip"
            echo ""
        fi
    fi
}

coreTestDir="core"
coreTests=(config_test.py lattice_test.py log_test.py utils_test.py planewave_test.py)
vaspTestDir="vasp"
vaspTests=(poscar_test.py)

for ct in ${coreTests[@]}
do
    fn=$coreTestDir/$ct
    run_one_unittest $fn
done

for vt in ${vaspTests[@]}
do
    fn=$vaspTestDir/$vt
    run_one_unittest $fn
done