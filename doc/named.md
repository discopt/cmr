# Named Matrices # {#named}

There exist several matrices that play a special role in theory.

The command

    cmr-named IN-MAT [OPTION]...

determines whether the [matrix](\ref file-formats-matrix) given in file `IN-MAT` is (up to row/column permutations) one of the named matrices below.

**Representation matrices of:**
  - `R_10`       regular matroid \f$ R_{10} \f$

**Variants (can be combined):**
  - `<NAME>*   ` refers to the dual matroid.
  - `<NAME>.<k>` refers to a specific representation matrix, \f$ k = 1,2, \dots \f$

**Other matrices:**
  - `I_<SIZE>`   Identity matrix of order SIZE.

**Options:**
  - `-i FORMAT`   Format of file `IN-MAT`; default: [dense](\ref dense-matrix).
  - `-e EPSILON`  Allows rounding of numbers up to tolerance EPSILON; default: \f$ 10^{-9} \f$.



Formats for matrices: [dense](\ref dense-matrix), [sparse](\ref sparse-matrix)

If `IN-MAT` is `-` then the [matrix](\ref file-formats-matrix) is read from stdin.

### Algorithm ###

The implemented algorithm is based on Section 18 of [A decomposition theory for matroids. V. Testing of matrix total unimodularity](https://doi.org/10.1016/0095-8956(90)90030-4) by Klaus Truemper (Journal of Combinatorial Theory, Series B, 1990).
For a matrix \f$ M \in \{-1,0,1\}^{m \times n} \f$ it runs in \f$ \mathcal{O}( \min(m^2 \cdot n, m \cdot n^2) ) \f$ time.

### C Interface ###

The corresponding function in the library is

  - CMRcamionComputeSigns() Computes a Camion-signed version of a given ternary matrix.

and is defined in \ref camion.h.
