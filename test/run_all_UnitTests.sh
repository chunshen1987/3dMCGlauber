#!/usr/bin/env bash

cp -r ../eps09 ./

for iexe in `ls | grep "unittest"`
do
    echo "$(tput bold)$(tput setaf 5)Running Unit Tests for $iexe ... $(tput setaf 0)"
    ./$iexe
done

rm -fr eps09
