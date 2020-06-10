
all: svd main

svd: 
	gcc -c src/svd.c

main: svd
	gcc -o a.out src/main.c svd.o -lm