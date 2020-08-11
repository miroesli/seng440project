# Custom Lookup Tables

## Requirements

- python3.8
- pip packages in [requirements.txt](../requirements.txt)

## Install

Use `pip install -r requirements.txt` in root to install requried packages.

## Configuration

Edit the global variables in the [create_lookup_table.py](../scripts/create_lookup_table.py) script according to your needs.

- `DEBUG` enables printing the result to the command line.
- `WRITE_TO_FILE` automatically updates the corresponding lookup table in [tables](../src/tables).
-
- `ARCTAN_TABLE_SIZE` defines the size of the arctangent lookup table.
- `SINCOS_TABLE_RANGE` defines the size of the sin and cos lookup table.
- `VALUES` defines the number of values to generate between 0 and the size of the lookup table.
- `SCALE_FACTOR` defines the scaling to apply to the float values for fixed-point notation.
- `SCALE_FACTOR_ARCTAN` the scale factor used for arctangent.
- `SCALE_FACTOR_SINCOS` the scale factor used for sine and cosine.

## Execution

To run the script run using python and specify one of `sin`, `cos`, or `arctan`.

```bash
python create_lookup_table.py arctan
```

You can also pipe it into the file of your choice if you have the `DEBUG` constant set to `True`.

```bash
python create_lookup_table.py arctan >> config.h
```
