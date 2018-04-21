SphericalmovingGRP.exe: SphericalmovingGRP.o VIPLimiter.o
	g++ -o SphericalmovingGRP.exe SphericalmovingGRP.o VIPLimiter.o -lm
SphericalmovingGRP.o: SphericalmovingGRP.cpp
	g++ -c SphericalmovingGRP.cpp
VIPLimiter.o: VIPLimiter.cpp
	g++ -c VIPLimiter.cpp

clean:
	rm -f SphericalmovingGRP.exe SphericalmovingGRP.o VIPLimiter.o
.PHONY:all clean
