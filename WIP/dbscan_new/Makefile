app: objects/main.o objects/dbscan.o objects/Points.o objects/InputReader.o objects/utils.o
	g++ objects/main.o objects/dbscan.o objects/Points.o objects/InputReader.o objects/utils.o -o app
objects/main.o: main.cpp
	g++ -c main.cpp -o objects/main.o
objects/dbscan.o: dbscan.cpp
	g++ -c dbscan.cpp -o objects/dbscan.o
objects/Points.o: Points.cpp
	g++ -c Points.cpp -o objects/Points.o
objects/InputReader.o: InputReader.cpp
	g++ -c InputReader.cpp -o objects/InputReader.o
objects/utils.o: utils.cpp
	g++ -c utils.cpp -o objects/utils.o
debug: main.cpp dbscan.cpp Points.cpp InputReader.cpp utils.cpp
	g++ -g main.cpp dbscan.cpp Points.cpp InputReader.cpp utils.cpp -o app
