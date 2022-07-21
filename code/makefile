ifndef CC
	CC = gcc
endif
ifndef CXX
	CXX = g++
else
	CXX += -std=c++17
endif
RM=rm -f
RMDIR=rm -rf

TRIMESH=dependency/trimesh2
RSTREE=dependency/rstree
EIGEN=dependency/eigen

INCS=-Iinclude/basic/ -Iinclude/c2o/ -Iinclude/4pcs/ -I$(RSTREE)/ -I$(TRIMESH)/include/ -I$(EIGEN)/
OBJS=$(RSTREE)/*.o
LIBS=-L$(TRIMESH)/lib.Linux64/ -ltrimesh -lgluit -fopenmp

## if cygwin is used, there are some missing functions causing the compilation to fail
## using below flag to solve the issue
# OPTS=-D_CYGWIN_WORKAROUND

SRC=src
OUT=out

.PHONY: clean dep_trimesh dep_rstree run_index run_query

all: run_index run_query

run_index: dep_trimesh dep_rstree mkdir_out
	$(CXX) $(INCS) $(OBJS) $(SRC)/RunIndex.cpp -o $(OUT)/run_index.out $(LIBS) $(OPTS)

run_query: dep_trimesh dep_rstree mkdir_out
	$(CXX) $(INCS) $(OBJS) $(SRC)/RunQuery.cpp -o $(OUT)/run_query.out $(LIBS) $(OPTS)

dep_trimesh:
	make -C $(TRIMESH)

dep_rstree:
	make -C $(RSTREE)

mkdir_out:
	mkdir -p $(OUT)/

clean:
	make clean -C $(TRIMESH)
	make clean -C $(RSTREE)
	$(RM) $(OUT)/*.out
