# Efficient Spatial Index for 3D Object Retrieval

This repository contains the code and the technical report for the submission "Efficient Spatial Index for 3D Object Retrieval".

## Requirements

* Linux environment
* GCC v10 (or another GCC version supporting C++17)

## Usage

1. Download this repository.

2. Change the working direction into code/.

3. Compile with "make" command.
   - Note that if you are using Cygwin, you need to uncomment line 22 in file "code/makefile" and uncomment line 39 in file "code/dependency/trimesh2/Makedefs.Linux".
   - After successful compilation, two executables will be generated under code/out/: run_index.out and run_query.out.
   - We also provide the two pre-generated Linux executables under "code/out/pre_compiled/" to be readily executed. To run them, just use "./run_pre_compiled.sh".

4. Run our test programs with "./run_all_tests.sh".
   - This will run C2O index building, C2O index loading and C2O querying for each of the 3 testing datasets:
     * testdata/object1/: 7 database objects from our dataset object and one of them as the query object.
     * testdata/object2/: 1000 database objects from our dataset object and a pre-generated random query point cloud in our experiments.
     * testdata/indoor/: 3 database objects from our dataset indoor and a pre-generated random query point cloud in our experiments.
   - As can be seen in run_all_tests.sh, the program basically take a dataset folder as the input parameter, e.g.,:

```sh
./out/run_index.out testdata/object1/
```

## Dataset Format

For each dataset, the detailed file formats are described below.

1. meta.txt: the metadata of a list of point clouds, where the first line is the number of the point clouds, and for each rest line the ID and filename (in this dataset folder) is shown. An example of 3 point clouds is shown below:
```
3
0 recon_bedroom.ply
1 recon_boardroom.ply
2 recon_loft.ply
```

2. Database point clouds: the database point clouds, which are basically a list of point clouds, can be specified by a meta.txt file. If meta.txt does not exist in a dataset folder, the database will contain all the point cloud files ending up with extension name "ply" (a popular 3D object file format, and the detailed description of a "ply" file can be found in http://paulbourke.net/dataformats/ply/)

3. Query point cloud: a point cloud file with extension name "ply". By default, this will be "query.ply".

## Configuration File Format

1. config_index.txt: configuration file for building index (each line corresponds to a parameter). The following list the possible parameters that can be specified:
   * r: the parameter r for building C2O index
   * index_path: the folder to place the generated index, default is the current dataset folder
   * index_filename: the filename of the generated index, default is "index"
   * show_progress_bar: 1 means that the progress bar will be shown, and 0 otherwise

   An example is shown below:
```
r=300
index_path=<default>
index_filename=<default>
show_progress_bar=1
```

2. config_query.txt: configuration file for executing query (each line corresponds to a parameter). The following list the possible parameters that can be specified:
   * delta: the delta value for querying (one unit is defined to be the diameter of the query point cloud)
   * query_path: the folder where the query point cloud locates, default is the current dataset folder
   * query_filename: the filename of the query point cloud, default is "query.ply"
   * index_path: the folder where the index file locates, default is the current dataset folder
   * index_filename: the filename of the index to be loaded, default is "index"
   
   An example is shown below:
```
delta=0.01
query_path=<default>
query_filename=00021.ply
index_path=<default>
index_filename=<default>
```

## Output Format

1. The program run_index.out outputs two files: index.ent and index.rst (both of them are in binary form).
   * index.ent: stores the information of index entries, i.e., the dababase ID and the four point IDs of each indexed relative-distance representation.
   * index.rst: stores the R*-tree.

2. The program run_query.out outputs the query result file: output.txt. In this file, the first line records the query execution time in second. The second line is the number of matched database objects. Then, for each matched object, the file shows the database ID, filename, the distance (one unit is defined to be the diameter of the query point cloud) from the query point cloud to the database matched object and the optimal transformation of the matching sequentially. An example is shown below:
```
0.016
2
985 08279s01.ply 0.00195051
-0.029842 0.605524 -0.795267 1656.66
-0.950676 0.228561 0.209702 2353.33
0.308746 0.762299 0.568836 1167.3
0 0 0 1
984 08279s00.ply 0.0019527
-0.029842 0.605524 -0.795267 1656.66
-0.950676 0.228561 0.209702 2353.33
0.308746 0.762299 0.568836 1167.3
0 0 0 1
```

   This example shows a query execution which takes 0.016s. Two database objects are returned, with their IDs 985 and 984 respectively, their filenames 08279s01.ply and 08279s00.ply respectively and their distances 0.00195051 and 0.0019527 respectively. After each ID and filename, the optimal transformation related to the matched object is shown in 4x4 matrix form.

## Additional Description

Additionally, we describe each folder in our code.

   * dependency/: the dependencies for this project. All the third-party code (i.e., trimesh2 and eigen) have been pre-included and integrated into the compilation process.
   * include/: the header files for using the C2O indexing and querying.
   * out/: the generated executable files.
   * src/: the source code of our C2O indexing and querying programs.
   * testdata/: the test datasets.
