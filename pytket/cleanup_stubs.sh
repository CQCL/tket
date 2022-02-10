#!/bin/bash

for file in $(ls -d pytket/_tket/*.pyi)
do
    echo "cleaning up $file"
    echo "from .type_aliases import *" > tmp_pyi
    cat $file >> tmp_pyi
    sed 's/numpy\.ndarray\[[^]]*\]\]/numpy.ndarray/g' tmp_pyi | grep -v -F -f stub_blacklist.txt > $file
done

cp binders/type_aliases.py pytket/_tket/type_aliases.pyi
rm tmp_pyi