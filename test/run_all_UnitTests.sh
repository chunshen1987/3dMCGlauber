#!/usr/bin/env bash

for iexe in `ls | grep "unittest"`
do
    echo "$(tput bold)$(tput setaf 5)Running Unit Tests for $iexe ... $(tput setaf 0)"
    ./$iexe
done
