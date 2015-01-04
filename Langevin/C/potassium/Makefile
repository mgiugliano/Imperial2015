CXX      = g++
CFLAGS   = -O2 -ffloat-store -Wall
LIBS     = -lm
SRC      = src
OBJECT_FILES = microscopic.o effective.o fox.o main.o

EXEC_FILES = potassium

all : $(EXEC_FILES)

rebuild : clean all

$(EXEC_FILES) : $(OBJECT_FILES)
	$(CXX) -o $(EXEC_FILES) $(OBJECT_FILES) $(LIBS)

main.o: $(SRC)/main.cpp $(SRC)/randlib.h
	$(CXX) $(CFLAGS) -c $(SRC)/main.cpp 

microscopic.o: $(SRC)/microscopic.cpp $(SRC)/microscopic.h $(SRC)/randlib.h
	$(CXX) $(CFLAGS) -c $(SRC)/microscopic.cpp 

effective.o: $(SRC)/effective.cpp $(SRC)/effective.h $(SRC)/randlib.h
	$(CXX) $(CFLAGS) -c $(SRC)/effective.cpp 

fox.o: $(SRC)/fox.cpp $(SRC)/fox.h $(SRC)/randlib.h
	$(CXX) $(CFLAGS) -c $(SRC)/fox.cpp 

micro : mainMicro.o effectiveMicro.o
	$(CXX) -o $(EXEC_FILES)Fast mainFast.o effectiveFast.o $(LIBS)

clean:
	rm -f $(SRC)/*~
	rm -f *~
	rm -f *.o

superclean: clean
	rm -f $(EXEC_FILES)

