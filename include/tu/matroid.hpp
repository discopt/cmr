#pragma once

#include <tu/config.h>
#include <tu/export.h>
#include <tu/matrix.hpp>

namespace tu
{ 

  /**
    * \brief Performs a binary pivot at \p row, \p column.
    *
    * Performs a binary pivot at \p row, \p column. This means
    * \f[ M_{i,j} = M_{i,j} - \frac{ M_{i, \text{\texttt{column}}} \cdot
    * M_{\text{\texttt{row}}, j} }{ M_{\text{\texttt{row}}, \text{\texttt{column}}} } \f].
    *
    * The implementation only ensures correctness for binary matrices. It flips all entries at i,j
    * for which i is different from \p row, j is different from \p column, and for which the
    * entries at i, \p column and \p row, j are nonzero. Flipping means to subtract the old value
    * from 1.
    */

  template <typename V>
  void matrixBinaryPivot(DenseMatrix<V>& matrix, Index row, Index column)
  {
    assert(matrix.get(row, column) != 0);

    for (Index i = 0; i < matrix.numRows(); ++i)
    {
      if (i == row)
        continue;

      const Value& first = matrix.get(i, column);
      if (first == 0)
        continue;

      for (Index j = 0; j < matrix.numColumns(); ++j)
      {
        if (j == column)
          continue;

        const Value& second = matrix.get(row, j);
        if (second == 0)
          continue;

        matrix.set(i, j, 1 - matrix.get(i, j));
      }
    }
  }

  /**
    * \brief Performs a ternary pivot at \p row, \p column.
    *
    * Performs a ternary pivot at \p row, \p column. This means
    * \f[ M_{i,j} = M_{i,j} - \frac{ M_{i, \text{\texttt{column}}} \cdot
    * M_{\text{\texttt{row}}, j} }{ M_{\text{\texttt{row}}, \text{\texttt{column}}} } \f].
    *
    * The implementation only ensures correctness for ternary matrices. It computes the resulting
    * values only based on the signs of the entries.
    */

  template <typename V>
  void matrixTernaryPivot(DenseMatrix<V>& matrix, Index row, Index column)
  {
    const Value& base = matrix.get(row, column);
    assert(base != 0);
    signed char denominator = base > 0 ? 1 : -1;

    for (Index i = 0; i < matrix.numRows(); ++i)
    {
      if (i == row)
        continue;

      const Value& first = matrix.get(i, column);
      if (first == 0)
        continue;
      signed char firstInt = first > 0 ? 1 : -1;

      for (Index j = 0; j < matrix.numColumns(); ++j)
      {
        if (j == column)
          continue;

        const Value& second = matrix.get(row, j);
        if (second == 0)
          continue;
        signed char numerator = firstInt * (second > 0 ? 1 : -1);

        const Value& value = matrix.get(i, j);
        signed char result = (int(get(i,j)) - numerator / denominator + 4) % 3 - 1;
        assert(result >= -1 && result <= 1);
        matrix.set(i, j,  result);
      }
    }
  }

  /**
   * \brief Matroid elements of a \ref Matrix.
   * 
   * Matroid elements of a \ref Matrix. Initially, rows will be labeled from -1 to -m and columns
   * from 1 to n, where the matrix is m-by-n.
   */

  class Matroid
  {
  public:
    /// Type for row/column labels.
    typedef short Element;

    /**
     * \brief Constructs matroid for a matrix.
     * 
     * Constructs matroid for an m-by-n matrix. Rows will be labeled from -1 to -m and columns will
     * be labeled from 1 to n.
     */

    template <typename M>
    Matroid(const M& matrix)
      : _rows(matrix.numRows()), _columns(matrix.numColumns())
    {
      for (short r = 0; r < _rows.size(); ++r)
        _rows[r] = -1 - r;
      for (short c = 0; c < _columns.size(): ++c)
        _columns[c] = 1 + c;
    }

    /**
     * \brief Returns the element corresponding to \p row.
     */

    Element row(std::size_t row) const
    {
      return _rows[row];
    }

    /**
     * \brief Returns the element corresponding to \p column.
     */

    Element column(std::size_t column) const
    {
      return _columns[column];
    }

  private:
    /// Row labels.
    std::vector<Element> _rows;
    /// Column labels.
    std::vector<Element> _columns;
  };

} /* namespace tu */
