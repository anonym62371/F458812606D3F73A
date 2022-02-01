#!/bin/bash

./out/pre_compiled/run_index.out testdata/object1/
./out/pre_compiled/run_index.out testdata/object2/
./out/pre_compiled/run_index.out testdata/indoor/

./out/pre_compiled/run_query.out testdata/object1/
./out/pre_compiled/run_query.out testdata/object2/
./out/pre_compiled/run_query.out testdata/indoor/