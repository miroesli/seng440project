# Compile Commands

## Basics

### Compile into assembly

_If arm-linux-gcc is not recognized, use `gcc` instead_

```bash
arm-linux-gcc -static -S file.c
```

add `-O3` for compiler optimization which does automated ufnction inlining

```bash
arm-linux-gcc -static -S -O3 file.c
```

To assemble modifications

```bash
arm-linux-gcc -static file.s
```

```bash
arm-none-linux-gnueabi-gcc -static -march=armv5 file.c -o file.exe
```

Compile floating point assembly

```bash
arm-linux-gcc -mfloat-abi=soft-S file.c
arm-linux-gcc -mfloat-abi=softfp-S file.c
```

### CPU profiling and Profile-driven compilation

Generate the executable file with profiling turned

```bash
gcc -pg file.c
```

Run the executable to produce the profile information (gmon.out file)

```bash
./a.out
```

Read the profile information

```bash
gprof
```

Analyse the profile informationa nd update the source file `file.c`

Recompile the source file with profiling turned off

```bash
gcc file.c
```

Another performance profiling is to use `perf`

```bash
gcc file.c -o file.exe
perf record ./file.exe
perf report
```

Use valgrind for various other profiling stats

```bash
gcc file.c -o file.exe
valgrind -tool=callgrind ./file.exe
callgrind_annotate callgrind.out.PID | grep function_1

# And to determine cache misses

gcc file.c -o file.exe
valgrind -tool=cachegrind -branch-sim=yes ./file.exe
```

### Instatiating machine-level custom instructions

Add and inline `__asm__ (...)` assmbly code / function directly in c code and then run compiler

```bash
arm-none-linux-gnueabi-gcc -static -march=armv5 -S file.c
```

To assemble use

```bash
arm-none-linux-gnueabi-gcc-static -march=armv5 file.s -o file
```

Simualting the resulting executable

```bash
qemu-arm file.exe
```

### Floating point operations

To analyze fixed point arithmetic, viewing floating point operations may be helpful

To emulate floating point operations run with `-mfloat-abi=soft`

```bash
arm-linux-gcc -mfloat-abi=soft -S file.c
```

However, it can be more helpful to allow the machine level operations with floating point support

```bash
arm-linux-gcc -mfloat-abi=softfp -S file.c
```
