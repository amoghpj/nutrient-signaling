model: pydstool2cpp.py
	python pydstool2cpp.py

# Machine specific
spindrift: main.cpp model.cpp 
	g++ -std=c++11 main.cpp model.cpp -o spindrift_main.o 


honoria: main.cpp model.cpp 
	g++ -std=c++11 main.cpp model.cpp -o honoria_main.o 

pirbright: main.cpp model.cpp 
	g++ -std=c++11 main.cpp model.cpp -o pirbright_main.o 

wyatt: main.cpp model.cpp 
	g++ -std=c++11 main.cpp model.cpp -o wyatt_main.o 

simmons: main.cpp model.cpp 
	g++ -std=c++11 main.cpp model.cpp -o simmons_main.o #20190226.o 

parsloe: main.cpp model.cpp 
	g++ -std=c++11 main.cpp model.cpp -o parsloe_main.o 

mnemosyne: main.cpp model.cpp 
	g++ -std=c++11 main.cpp model.cpp -o mnemosyne_main.o 

cuthbert: main.cpp model.cpp 
	g++ -std=c++11 main.cpp model.cpp -o cuthbert_main.o 

carmody: main.cpp model.cpp 
	g++ -std=c++11 main.cpp model.cpp -o carmody_main.o 


all: main.cpp model.cpp #simulator.cpp
	g++ -std=c++11 main.cpp model.cpp -o main.o 

sim: main.o
	./main.o --ode

plot: main.o 
	./main.o --ode --plot

clean:
	rm *.dat *.o
