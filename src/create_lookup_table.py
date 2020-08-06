""" Script to create a lookup table array for C

The user provides an input between 0 and MAX_INPUT size in C
If outside bounds, use +-pi/2

Based on trig function supplied run accordingly
"""

import sys
import numpy as np

# Size of the lookup table -> Integer bit size 13, 14th bit is used for sign in C code
MAX_INPUT = 1 << 13
# Range of lookup table
LOOKUP_TABLE_SIZE = 4096  # MAX_INPUT
# scale factor between float and fixed point integer
SCALE_FACTOR = 1 << 30
# Precision of the lookup tables
PRECISION = 5

# TODO add precision to tan (between integers) and cos and sin (how many decimal places of the range 0-1)

"""Create lookup table for arcsin

Input range for the lookup table will be between 0 and 1
"""


def create_arcsin_lookup_table():
    pass


"""Create lookup table for arccos

Input range for the lookup table will be between 0 and 1
"""


def create_arccos_lookup_table():
    pass


"""Create lookup table for arctan

Input range is predefined based on precision desired.
Maximum range is between 0 and integer size.
"""

# Consider creating and approximator instead for tan, but for cos and sin just use lookup table.


def create_arctan_lookup_table():
    for index in range(LOOKUP_TABLE_SIZE):
        # use selected trig function
        val = np.arctan(index)*SCALE_FACTOR
        try:
            print("%d" % (val), end="")
        except Exception:
            print("NaN", end="")

        # print new lines
        if (index + 1) % 10 == 0:
            print("")
        else:
            print(" ", end="")


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

    switcher = {
        "sin": create_arcsin_lookup_table,
        "cos": create_arccos_lookup_table,
        "tan": create_arctan_lookup_table,
    }
    output_lookup_table = switcher.get(trig_function, np.arctan)

    # Create c code lookup definition
    print("static const __flash uint32_t arc%s_lookup[%d] = "
          % (trig_function, LOOKUP_TABLE_SIZE))
    print("{")
    output_lookup_table()
    print("};")


if __name__ == "__main__":
    main()
