/*
 *
 * main.c
 *
 */

#include <stdio.h>
volatile int c;

int add(int a, int b) { return a + b; }

int main(void)
{
    int a = 1, b = 2;
    c = add(a, b);
    printf("a + b = %i\n", c);
    exit(1);
}