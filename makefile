optimization = -O3 -g -Wno-deprecated

iHS_calc: PData.cc Ehh.cc
	g++ $(optimization) -o iHS_calc PData.cc Ehh.cc

clean:
	rm iHS_calc
