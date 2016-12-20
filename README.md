# graphcut
Graphcut with the Boykov Kolmogorov, with app made in GTK+ 3. It also uses Armadillo.

To compile:

g++ -std=c++14 `pkg-config --cflags gtk+-3.0` -o main main.cpp `pkg-config --libs gtk+-3.0` -larmadillo
