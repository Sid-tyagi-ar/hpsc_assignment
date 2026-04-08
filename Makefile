all:
	g++ -O3 -fopenmp HPSC_Assignment.cpp -o sim

run: all
	./sim

plot:
	python3 plot_helper.py

clean:
	rm -f sim *.txt *.png
