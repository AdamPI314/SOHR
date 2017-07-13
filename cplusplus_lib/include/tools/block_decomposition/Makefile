ICC= mpic++
OPT= -Ofast

P= 4

SRC= main.cpp block_decomposition.cpp
EXE= main

all:test

test: $(SRC)
	$(ICC) $(OPT) $(SRC) -o $(EXE) -lboost_mpi -lboost_serialization

run:
	mpirun -np $(P) ./$(EXE)

clean:
	rm $(EXE)
