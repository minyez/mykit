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
                exit $status
            fi
            echo ""
        else
            echo "!!WARN: $fn not found. Skip"
            echo ""
        fi
    fi
}

coreTestDir="core"
coreTests=(config_test.py
           bandstructure_test.py
           cell_test.py
           log_test.py
           ion_test.py
           utils_test.py
           planewave_test.py
           kmesh_test.py
           xc_test.py
           control_test.py
           symmetry_test.py)

vaspTestDir="vasp"
vaspTests=(poscar_test.py
           incar_test.py
           kpoints_test.py
           potcar_test.py
           xml_test.py)

w2kTestDir="wien2k"
w2kTests=(inputs_test.py struct_test.py)

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

for wt in ${w2kTests[@]}
do
    fn=$w2kTestDir/$wt
    run_one_unittest $fn
done
