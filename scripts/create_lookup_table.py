""" Script to create a lookup table array for C

A table for sin, cos or arctan is generated based on selection.
"""

import sys
import numpy as np

# DEBUG
DEBUG = False
# Print to file
WRITE_TO_FILE = True
# Range of arctan lookup table
ARCTAN_TABLE_RANGE = 10
# Constant range
SINCOS_TABLE_RANGE = np.pi
# the resolution of the table
VALUES = 10000
# scale factor between float and fixed point integer
SCALE_FACTOR_ARCTAN = 1 << 13
SCALE_FACTOR_SINCOS = 1 << 13

# TODO Consider creating and approximator instead for tan,
# but for cos and sin just use lookup table.

"""Create lookup table
arctan - call lookup with value [0, MAX_VALUE]
Note: requries updating sign after result.
sin/cos - call lookup with value [0, 2pi]
Note: must modulo before providing input

When accessing the table the user must multiply by values and divide
by range to scale for table size and resolution.

This function prints the results in 10 columns to be
stored in a c header file.
"""


def create_lookup_table(f, value_function, range, values):
    for index, x in enumerate(np.arange(0, range, range/values)):
        y = value_function(x)
        # TODO values not being rounded correctly - this is fine for fixed point?
        if DEBUG:
            print('{:= 11d}'.format(int(y)), end="")
        if WRITE_TO_FILE:
            f.write('{:= 11d}'.format(int(y)))
        if (index + 1) % 10 == 0:
            if DEBUG:
                print(",")
            if WRITE_TO_FILE:
                f.write(",\n")
        else:
            if DEBUG:
                print(", ", end="")
            if WRITE_TO_FILE:
                f.write(", ")


"""Print the correct usage of the script."""


def usage():
    print("Please provide the trig function you'd like to create a lookup table for.")
    print("eg: python create_lookup_table.py [sin|cos|arctan]")


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
        "arctan": lambda y: np.arctan(y) * SCALE_FACTOR_ARCTAN,
        "sin": lambda y: np.sin(y) * SCALE_FACTOR_SINCOS,
        "cos": lambda y: np.cos(y) * SCALE_FACTOR_SINCOS
    }
    value_function = value_func_switcher.get(trig_function)

    if trig_function == "arctan":
        range = ARCTAN_TABLE_RANGE
    else:
        range = SINCOS_TABLE_RANGE

    # Create c code lookup definition
    f = None
    if WRITE_TO_FILE:
        f = open("../src/tables/%s_lookup_table.h" % trig_function, "w")

    header_content = \
        """/*
 *
 * %s_lookup_table.h
 *
 */

#ifndef %s_lookup_table_h
#define %s_lookup_table_h

#include "../svd_math.h"

"""
    if DEBUG:
        print(header_content % (trig_function, trig_function, trig_function))
    if WRITE_TO_FILE:
        f.write(header_content % (trig_function, trig_function, trig_function))
    if DEBUG:
        print("static const fixed_point_t %s_lookup_table[%d] = "
              % (trig_function, VALUES))
    if WRITE_TO_FILE:
        f.write("static const fixed_point_t %s_lookup_table[%d] = \n"
                % (trig_function, VALUES))
    if DEBUG:
        print("{")
    if WRITE_TO_FILE:
        f.write("{\n")

    # write the numbers
    create_lookup_table(f, value_function, range, VALUES)

    if DEBUG:
        print("}; \n\n#endif")
    if WRITE_TO_FILE:
        f.write("}; \n\n#endif")

    if WRITE_TO_FILE:
        print("Wrote lookup table to %s_lookup_table.h" % trig_function)
    if WRITE_TO_FILE:
        f.close()


if __name__ == "__main__":
    main()
