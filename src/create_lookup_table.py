""" Script to create a lookup table array for C

The user provides an input between 0 and MAX_INPUT size in C
If outside bounds, use +-pi/2

Based on trig function supplied run accordingly
"""

import sys
import numpy as np

# Range of arctan lookup table
ARCTAN_TABLE_RANGE = 30  # MAX_INPUT
SINCOS_TABLE_RANGE = 2 * np.pi
# the resolution of the table
VALUES = 300
# scale factor between float and fixed point integer
SCALE_FACTOR = 1 << 31

# TODO Consider creating and approximator instead for tan,
# but for cos and sin just use lookup table.

"""Create lookup table
arctan - call lookup with value [0, MAX_VALUE]
(requries updating sign after result)

must multiply by values and divide by range to scale for table

---

input range of arccos is [-1,1] but add 1 to scale to [0,2]
when accessing lookup table. Input range for the lookup table
will be between 0 and 2 scaled to an integer based on precision

Input range is predefined based on precision desired.
Maximum range is between 0 and integer size.

If higher precision is used, assuming input is scaled accordingly
to be an integer.
"""


def create_lookup_table(value_function, range, values):
    for index, x in enumerate(np.arange(0, range, range/values)):
        # print(x, index)
        y = value_function(x)
        # print the value
        # print("%010d" % (y), end="")
        print('{:= 11d}'.format((int(y))), end="")
        if (index + 1) % 10 == 0:
            print(",")
        else:
            print(", ", end="")


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
        if trig_function not in ["sin", "cos", "arctan"]:
            usage()
            exit(1)

    value_func_switcher = {
        # * SCALE_FACTOR
        "arctan": lambda y: np.arctan(y) * SCALE_FACTOR,
        "sin": lambda y: np.sin(y) * SCALE_FACTOR,
        "cos": lambda y: np.cos(y) * SCALE_FACTOR
    }
    value_function = value_func_switcher.get(trig_function)

    if trig_function == "arctan":
        range = ARCTAN_TABLE_RANGE
    else:
        range = SINCOS_TABLE_RANGE

    # Create c code lookup definition
    print("static const fixed_point_t %s_lookup[%d] = "
          % (trig_function, VALUES))
    print("{")
    create_lookup_table(value_function, range, VALUES)
    print("};")


if __name__ == "__main__":
    main()
