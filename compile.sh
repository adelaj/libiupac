gcc -std=c11 -Wall -O3 -mavx -o ./libiupac.o -c ./Source/libiupac.c
gcc -std=c11 -Wall -O3 -mavx -o ./main.o -I./Source -c ./Source-Test/main.c
gcc -o ./libiupac-test ./main.o ./libiupac.o
