
all: svd svd_math main
test:  svd_math test_c

svd: 
	gcc -c src/svd.c
svd_math:
	gcc -c src/svd_math.c
main: svd svd_math
	gcc -o a.out src/main.c svd.o svd_math.o -lm
test_neon:
	gcc -mfpu=neon -O3 -S src/test_neon.c
	gcc -mfpu=neon -O3 -o file.exe src/test_neon.c
test_c: svd_math
	gcc -o test.out src/test.c svd_math.o -lm
