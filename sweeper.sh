#!/usr/bin/env bash

results_folder=$1

mkdir -p $results_folder

mv string* $results_folder
mv events_summary.dat $results_folder
cp input $results_folder
./utilities/combine_events_into_hdf5.py $results_folder
