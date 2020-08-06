# Custom Lookup Tables

## Requirements

- python3.8
- pip packages in [requirements.txt](../requirements.txt)

## Install

Use `pip install -r requirements.txt` in root to install requried packages.

## Configuration

Edit the global variables in the [create_lookup_table.py](../src/create_lookup_table.py) script according to your needs.

- `ARCTAN_TABLE_SIZE` defines the size of the arctan lookup table.
- `SCALE_FACTOR` defines the scaling to apply to the float values before returning as integers.
- `PRECISION` defines the precision of the values. Maximum is 13 bits of info so

## Execution

To run the script run using python and specify one of `sin`, `cos`, or `tan`.

```bash
python create_lookup_table.py tan
```

You can also pipe it into the file of your choice

```bash
python create_lookup_table.py tan >> config.h
```
