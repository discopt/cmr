# Change Log # {#changes}

  - Revised the matroid decomposition structure.
  - Fixed a bug in the code for printing dot-files for transposed graphic/network matrices; reported by Christopher Hojny.
  - Added a `timeLimit` parameter to all potentially time intensive functions.
  - Bugfix for the `timeLimit` parameter.
  - Added \ref CMRdecIsSeriesParallelReduction and \ref CMRdecIsUnknown to interface for regular matroid decompositions.

## Version 1.3 ##

  - Redesigned as [Combinatorial Matrix Recognition](\ref index) Library.
  - Added recognition of [(co)graphic](\ref graphic) and [(co)network](\ref network) matrices.
  - Added recognition of [series-parallel matrices](\ref series-parallel).
  - Moved all C++ code to internal implementation. The library interface is only in C.
  - Added generators for network, graphic and random matrices.

## Version 1.2 ##

### 1.2h (2019-09-27) ###

  - Bugfix in search for violating submatrix.
  - Added direct test for regularity of binary matroids.

### 1.2g (2019-08-15) ###

  - Added support for complement total unimodularity.

### 1.2f (2017-03-14) ###

  - Raise an error message if an integer overflow occurs during k-modularity check; reported by Filippo Quondam.

### 1.2e (2016-02-28) ###

  - Bugfix for violator search; reported by Bala Krishnamoorthy.

### 1.2d (2015-02-06) ###

  - Bugfix for memory bug; reported by Tobias Windisch.

### 1.2c (2014-05-08) ###

  - Compatibility fix for more recent boost library; fix by Volker Braun.

### 1.2b (2012-10-12) ###

  - Compatibility fix for Mac OS; reported by Yvon Sauvagneau.

### 1.2 (2012-04-13) ###

  - Official version for article [Implementation of a unimodularity test](https://doi.org/10.1007/s12532-012-0048-x) by Matthias Walter and Klaus Truemper (Mathematical Programming Computation, 2013).
  - Polymake users can now analyze the matroid decomposition.

## Version 1.1 (2012-02-09) ##

  - Bug in graphicness test. Older versions may produce wrong results; reported by Paulo Seijas.
  - Update for polymake version 2.11 by Marc Pfetsch.
  
## Version 1.0 (2011-06-13) ##

  - First release of the software.

\anchor TUtest
# TUtest Library #

This library is the successor of the **TUtest** library, originally developed by Matthias Walter and Klaus Truemper for the recognition of [totally unimodular matrices](\ref tu).
The new name and additional functionality are available since version 1.3.
