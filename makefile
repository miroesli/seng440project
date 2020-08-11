
all: svd svd_math main
test:  svd_math test_c

svd: 
	gcc -mfpu=neon -O3 -c src/svd.c
svd_math:
	gcc -c src/svd_math.c
main: svd svd_math
	gcc -o a.out src/main.c svd.o svd_math.o -lm
test_neon:
	gcc -mfpu=neon -O3 -S src/tests/test_neon.c
	gcc -mfpu=neon -O3 -o file.exe src/tests/test_neon.c
test_c: svd_math
	gcc -mfpu=neon -O3 -o test.out src/tests/test.c svd_math.o -lm
