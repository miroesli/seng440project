
all: svd svd_math main
test: svd svd_math test_c

svd: 
	gcc -c src/svd.c
svd_math:
	gcc -c src/svd_math.c

main: svd svd_math
	gcc -o a.out src/main.c svd.o svd_math.o -lm
test_c: svd svd_math
	gcc -o test src/test.c svd.o svd_math.o -lm