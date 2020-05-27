/*
 *
 * main.c
 *
 */

#include "config.h"
#include "svd.h"

#include <stdio.h>
#include <stdlib.h>

volatile int c;

/* Template code */
int add(int a, int b) { return a + b; }

int main(void)
{
    int a = 1, b = 2;
    c = add(a, b);
    printf("a + b = %i\n", c);
    exit(1);
}