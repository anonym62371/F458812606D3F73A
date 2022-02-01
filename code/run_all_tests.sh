#!/bin/bash

./out/run_index.out testdata/object1/
./out/run_index.out testdata/object2/
./out/run_index.out testdata/indoor/

./out/run_query.out testdata/object1/
./out/run_query.out testdata/object2/
./out/run_query.out testdata/indoor/