app: kmeans.o kmeans_utils.o
	g++ kmeans.o kmeans_utils.o -o app
kmeans.o: kmeans.cpp
	g++ kmeans.cpp -c
kmeans_utils.o: kmeans_utils.cpp
	g++ kmeans_utils.cpp -lm -c