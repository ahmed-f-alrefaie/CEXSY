find . -name "*.cpp" | xargs icc -c -ipo -O3 -no-prec-div -std=c++0x -xHost; ls *.o | xargs icc -o main_vec.x -ipo -O3 -xHost -no-prec-div -std=c++0x -lpthread
