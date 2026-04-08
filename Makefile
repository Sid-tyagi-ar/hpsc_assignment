all:
	g++ -O3 -fopenmp HPSC_Assignment.cpp -o sim

run: all
	./sim

plot:
	python3 plot.py
	python3 combiner.py
clean:
	rm -f sim *.txt *.png
