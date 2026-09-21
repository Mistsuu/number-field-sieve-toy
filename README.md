# Number Field Sieve Toy Implementation

## Hello!

A very simple implementation that aids the understanding of the **Number Field Sieve Algorithm** used to factor numbers :-) Implemented from this paper: [link](https://personal.math.vt.edu/brown/doc/briggs_gnfs_thesis.pdf).

This **DOESN'T INCLUDE** optimizations used in the general programs implementing this algorithm, such as:
- Avoid computing square root of an element in `Z[O]` because the coeffients then will be too big.
- The Lanczos Algorithm used to compute solutions to the *matrix - vector equation in `GF(2)`*: `Ax = 0` in optimal space.
- *i don't really remember anything else since i just remember the main stuffs from the paper i'm currently reading... :'3*

Anyways, that means that this implementation is really slow, and really space consuming for some big parameters configurations...

## Prerequisites
- `python`
- `sage`

## Usage
Run
```python
# debug logs
DEBUG=1 python3 main.py

# no logs
python3 main.py
```


