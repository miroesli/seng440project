# Seng440 project

A project on Singular Value Decomposition (SVD) and optimizing its performance on an arm machine.

## Requirements

- qemu
- qemu-kvm
- qemu-system-arm
- libvirt-clients
- libvirt-daemon-system
- bridge-utils
- virt-manager
  - arch: arm
  - machine type: virt-2.11
  - Fedora 29
  - 2GB Ram
  - 1 CPU
- gcc compiler

## Installation

Create a VM in virt-manager using the settings specified in the [requirements](#requirements) section.

Add the following kernel arguments:

```bash
console=ttyAMA0
rw
root=LABEL=_/
rootwait
ipv6.disable=1
```

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

arm-none-linux-gnueabi-gcc -static -march=armv5 file.c -o file.exe

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

## Authors

Robert Tulip
Michail Roesli
