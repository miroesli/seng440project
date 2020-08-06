""" Script to create a lookup table array for C

The user provides an input between 0 and MAX_INPUT size in C
If outside bounds, use +-pi/2

Based on trig function supplied run accordingly
"""

import sys
import numpy as np
import math

# Size of the lookup table -> Integer bit size 13, 14th bit is used for sign in C code
MAX_INPUT = 1 << 13
# Range of arctan lookup table
ARCTAN_TABLE_SIZE = 4096  # MAX_INPUT
# scale factor between float and fixed point integer
SCALE_FACTOR = 1 << 30
# Precision of the lookup tables
PRECISION = 3

"""Print value and space or newline based on index in lookup table"""


def print_value(val, index):
    print("%d" % (val), end="")
    if (index + 1) % 10 == 0:
        print("")
    else:
        print(" ", end="")


# TODO Consider creating and approximator instead for tan,
# but for cos and sin just use lookup table. Also add precision
# to tan (between integers).

"""Create lookup table
input range of arccos is [-1,1] but add 1 to scale to [0,2]
when accessing lookup table. Input range for the lookup table
will be between 0 and 2 scaled to an integer based on precision

Input range is predefined based on precision desired.
Maximum range is between 0 and integer size.

If higher precision is used, assuming input is scaled accordingly
to be an integer.
"""


def create_lookup_table(value_function, size):
    for index in range(size):
        val = value_function(index)
        print_value(val, index)


"""Print the correct usage of the script."""


def usage():
    print("Please provide the trig function you'd like to create a lookup table for.")
    print("eg: python create_lookup_table.py [sin|cos|tan]")


"""Main function to evaluate which lookup table to create"""


def main():
    if len(sys.argv) != 2:
        usage()
        exit(1)
    else:
        trig_function = sys.argv[1]
        if trig_function not in ["sin", "cos", "tan"]:
            usage()
            exit(1)

    value_func_switcher = {
        "sin": lambda value: np.arcsin((value / pow(10, PRECISION))-1)
        * SCALE_FACTOR,
        "cos": lambda value: np.arccos((value / pow(10, PRECISION))-1)
        * SCALE_FACTOR,
        "tan": lambda value: np.arctan(value)
        * SCALE_FACTOR,
    }
    value_function = value_func_switcher.get(trig_function)

    lookup_table_size = ARCTAN_TABLE_SIZE if trig_function == "tan" \
        else pow(10, PRECISION)*2

    # Create c code lookup definition
    print("static const __flash uint32_t arc%s_lookup[%d] = "
          % (trig_function, lookup_table_size))
    print("{")
    create_lookup_table(value_function, lookup_table_size)
    print("};")


if __name__ == "__main__":
    main()
