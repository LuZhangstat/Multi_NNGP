
all:
	g++ -c -fpic -fopenmp nn.cpp -o nn.o
	g++ -shared -o nn.so nn.o	

clean:
	rm *.o *.so
