find . -name "*.cpp" | xargs icc -c -O0 -g -std=c++0x -xHost; ls *.o | xargs icc -o main_debug.x -O0 -g -std=c++0x -lpthread
