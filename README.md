 # Invariant-generating algorithm – Matrix invariants for program equivalence in LCTRSs

This artifact contains two implementations (one in GNU Octave and one in Python/SymPy) of the invariant-generating algorithm, as presented in the WPTE 2023 workshop paper:

**Matrix invariants for program equivalence in LCTRSs -  Kasper Hagens and Cynthia Kop.**
(Handle: https://hdl.handle.net/2066/312790)


## Contents

The artifact contains two implementations of the invariant-generating algorithm:

* `Python/`: Python/SymPy implementation
  * `invariants.py`: implementation of the invariant-generating algorithm and example computation corresponding to Example 6 in the paper

* `Octave/`: GNU Octave implementation
  * `invariants.m`: implementation of the invariant-generating algorithm
  * `demo.m`: example computation corresponding to Example 6 in the paper
  * `substitutionSumtwothree.m`: divergence substitution used in the example


## Algorithm 

The algorithm in both the GNU Octave and Python implementation computes all polynomial invariants of degree at most `d` for a given divergence pattern, characterized by the initialization vector and divergence substitution.

### Input: 
The required inputs are:

* `n`: the number of divergence variables;
* `d`: the degree of the polynomial invariants;
* `chi`: the divergence substitution;
* `init`: the initialization vector, which should be of length `n`.

Note: the initialization variables that occur in the theoretical description in the workshop paper do not need to be declared in the implementation, because the relevant initialization variables are  derived from the required input-argument `init`.

### Output: 
The algorithm creates divergence variables `V_div = {Y1, ..., Yn}` and outputs a generating set of polynomials: for each polynomial `p` in this generating set, the equation `p = 0` defines a polynomial invariant of `V_div`-degree at most `d`.

Both implementations also display intermediate objects used in the computation, including the monomial vector, the divergence matrix, and a basis for its kernel.

## Python implementation

### System requirements

The Python implementation requires:

* Python 3
* SymPy

Python 3.14.6 and SymPy 1.14.0 were used for the computations in this artifact. Other versions may also work, but have not been tested.

### How to use the `invariants` function

The `invariants` function is called using:

```python
invariants(n, d, chi, init)
```

Here:

- `n` is a positive integer 
- `d` is a non-negative integer
- `chi` is a Python function
- `init` is an `n x 1` `SymPy Matrix` (i.e. a size-`n` symbolic column vector)

The divergence substitution `chi` should be defined using 1-based indexing (i.e. using `Y[1], ..., Y[n]` rather than `Y[0], ..., Y[n-1]`). The function `index_zero_shift` takes care of converting the 1-based indexing to the 0-based indexing used internally by Python.

### Example

The file `invariants.py` contains a Python implementation of the invariant-generating algorithm as well as the computation for Example 6 in the paper, which considers invariant generation for the equivalence `sum2(x) ≈ sum3(x)`.

The example computation is executed automatically when `invariants.py` is run as a script:

```bash
python3 invariants.py
```


## GNU Octave implementation

### System requirements

The GNU Octave implementation requires:

* GNU Octave
* The GNU Octave symbolic package

The [GNU Octave symbolic package](https://octave.sourceforge.io/symbolic/) can be downloaded from the GNU Octave website.

Version 3.0.1 of the symbolic package was used for the computations in this artifact. Other versions may also work, but have not been tested.

The symbolic package must be installed and loaded in GNU Octave.

### How to use the `invariants` function

The `invariants` function is called using the following in the GNU Octave command window:

```octave
invariants(n, d, @chi, init)
```

Here:

- `n` is a positive integer 
- `d` is a non-negative integer
- `@chi` is a function handle referring to the function `chi`, which should be defined in a file named `chi.m` and placed in the same folder as `invariants.m`
- `init` is a symbolic vector. Such a vector can be declared in GNU Octave using the `syms` and `sym` commands. For example,

    ```octave
    syms x y;
    init = sym([x^2-1, x*y]);
    ```

    first creates symbolic variables `x` and `y`, which can then be used to declare the symbolic vector `(x^2-1, x*y)`.

### Example

The file `demo.m` contains a GNU Octave script that executes the invariant computation for Example 6 in the paper, which considers invariant generation for the equivalence `sum2(x) ≈ sum3(x)`.

The script sets the input parameters for this example and calls the `invariants` function using the divergence substitution defined in `substitutionSumtwothree.m`. As explained above, the files `demo.m`, `invariants.m`, and `substitutionSumtwothree.m` should be placed in the same folder before executing `demo.m`.
