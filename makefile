
all: svd main

svd: 
	gcc -c src/svd.c
svd_math:
	gcc -c src/svd_math.c

main: svd svd_math
	gcc -o a.out src/main.c svd.o svd_math.o -lm