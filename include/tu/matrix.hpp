#pragma once

#include <boost/numeric/ublas/matrix.hpp>

#include "matrix_transposed.hpp"
#include <iomanip>

namespace tu
{
  template <typename Iterator>
  class Range
  {
  public:
    Range(Iterator begin, Iterator end)
      : _begin(begin), _end(end)
    {

    }

    Range(const Range<Iterator>& other)
      : _begin(other._begin), _end(other._end)
    {

    }

    Iterator begin()
    {
      return _begin;
    }

    Iterator end()
    {
      return _end;
    }

  private:
    Iterator _begin;
    Iterator _end;
  };

  /**
   * \brief Base class for matrices whose entries have type \p V.
   */

  template <typename V>
  class Matrix
  {
  public:
    typedef V Value; /// Type of entries
    typedef std::size_t Index; /// Row/column index type

    struct Nonzero
    {
      Index row;
      Index column;
      Value value;

      Nonzero(Index r, Index c, Value v)
        : row(r), column(c), value(v)
      {

      }
    };
  };

  template <typename V>
  class DenseMatrix : public Matrix<V>
  {
  public:
    typedef typename Matrix<V>::Value Value;
    typedef typename Matrix<V>::Index Index;
    typedef typename Matrix<V>::Nonzero Nonzero;

    template <bool Row>
    class NonzeroIterator
    {
    public:
      NonzeroIterator(const DenseMatrix<V>& matrix, Index major)
        : _matrix(matrix), _major(major), _minor(0)
      {
        if (Row)
        {
          while (_minor < _matrix.numColumns() && _matrix.get(_major, _minor) == 0)
            ++_minor;
        }
        else
        {
          while (_minor < _matrix.numRows() && _matrix.get(_minor, _major) == 0)
            ++_minor;
        }
      }

      NonzeroIterator(const DenseMatrix<V>& matrix, Index major, int dummy)
        : _matrix(matrix), _major(major)
      {
        if (Row)
          _minor = _matrix.numColumns();
        else
          _minor = _matrix.numRows();
      }

      NonzeroIterator<Row>& operator++()
      {
        ++_minor;
        if (Row)
        {
          while (_minor < _matrix.numColumns() && _matrix.get(_major, _minor) == 0)
            ++_minor;
        }
        else
        {
          while (_minor < _matrix.numRows() && _matrix.get(_minor, _major) == 0)
            ++_minor;
        }
      }

      Nonzero operator*() const
      {
        if (Row)
        {
          assert(_major < _matrix.numRows());
          assert(_minor < _matrix.numColumns());
          return Nonzero(_major, _minor, _matrix.get(_major, _minor));
        }
        else
        {
          assert(_minor < _matrix.numRows());
          assert(_major < _matrix.numColumns());
          return Nonzero(_minor, _major, _matrix.get(_minor, _major));
        }
      }

      bool operator!=(const NonzeroIterator<Row>& other) const
      {
        assert(&_matrix == &other._matrix);
        assert(_major == other._major);
        return _minor != other._minor;
      }

    private:
      const DenseMatrix<V>& _matrix;
      std::size_t _major, _minor;
    };

    typedef Range<NonzeroIterator<true>> NonzeroRowRange;
    typedef Range<NonzeroIterator<false>> NonzeroColumnRange;
    
    /**
     * \brief Default constructor for 0x0 matrix.
     */

    DenseMatrix()
      : _data(nullptr), _numRows(0), _numColumns(0)
    {

    }

    /**
     * \brief Move constructor that steals the data from \p other.
     */

    DenseMatrix(DenseMatrix<V>&& other)
      : _data(other._data), _numRows(other._numRows), _numColumns(other._numColumns)
    {
      other._data = nullptr;
      other._numRows = 0;
      other._numColumns = 0;
    }

    /**
     * \brief Copy constructor.
     */

    template <typename V2>
    DenseMatrix(const DenseMatrix<V2>& other)
      : _numRows(other._numRows), _numColumns(other._numColumns)
    {
      std::size_t size = other._numRows * other._numColumns;
      _data = new Value[size];
      for (std::size_t i = 0; i < size; ++i)
        _data[i] = other._data[i];
    }

    /**
     * \brief Move assignment operator that steals the data from \p other.
     */

    DenseMatrix<V>& operator=(DenseMatrix<V>&& other)
    {
      _data = other._data;
      _numRows = other._numRows;
      _numColumns = other._numColumns;
      other._data = nullptr;
      other._numRows = 0;
      other._numColumns = 0;
      return *this;
    }

    /**
     * \brief Assignment operator.
     */

    template <typename V2>
    DenseMatrix<V>& operator=(const DenseMatrix<V2>& other)
    {
      _numRows = other._numRows;
      _numColumns = other._numColumns;
      std::size_t size = other._numRows * other._numColumns;
      _data = new Value[size];
      for (std::size_t i = 0; i < size; ++i)
        _data[i] = other._data[i];
      return *this;
    }

    /**
     * \brief Destructor.
     */

    ~DenseMatrix()
    {
      if (_data != nullptr)
        delete[] _data;
    }

    /**
     * \brief Constructs zero matrix of given size.
     */

    DenseMatrix(Index numRows, Index numColumns)
      : _numRows(numRows), _numColumns(numColumns)
    {
      std::size_t size = numRows * numColumns;
      _data = new Value[size];
      for (std::size_t i = 0; i < size; ++i)
        _data[i] = 0;
    }

    /**
     * \brief Returns the number of rows.
     */

    Index numRows() const
    {
      return _numRows;
    }

    /**
     * \brief Returns the number of columns.
     */

    Index numColumns() const
    {
      return _numColumns;
    }

    /**
     * \brief Returns entry at \p row, \p column.
     */

    const Value& get(Index row, Index column) const
    {
      assert(row < _numRows);
      assert(column < _numColumns);
      return _data[_numColumns * row + column]; 
    }

    /**
     * \brief Sets entry at \p row, \p column to copy of \p value.
     */

    template <typename V2>
    void set(Index row, Index column, const V2& value)
    {
      assert(row < _numRows);
      assert(column < _numColumns);
      _data[_numColumns * row + column] = value;
    }

    /**
     * \brief Returns the number of nonzeros in \p row.
     */

    std::size_t countRowNonzeros(Index row) const
    {
      std::size_t begin = _numColumns * row;
      std::size_t end = _numColumns * (row + 1);
      std::size_t count = 0;
      for (std::size_t i = begin; i < end; ++i)
      {
        if (_data[i] != 0)
          ++count;
      }
      return count;
    }

    /**
     * \brief Returns the number of nonzeros in \p column.
     */

    std::size_t countColumnNonzeros(Index column) const
    {
      std::size_t begin = column;
      std::size_t end = _numRows * _numColumns;
      std::size_t count = 0;
      for (std::size_t i = begin; i < end; i += _numColumns)
      {
        if (_data[i] != 0)
          ++count;
      }
      return count;
    }

    /**
     * \brief Returns a range for iterating over the nonzeros of \p row.
     */

    NonzeroRowRange iterateRowNonzeros(Index row) const
    {
      return NonzeroRowRange(NonzeroIterator<true>(*this, row),
        NonzeroIterator<true>(*this, row, 0));
    }

    /**
     * \brief Returns a range for iterating over the nonzeros of \p column.
     */

    NonzeroColumnRange iterateColumnNonzeros(Index column) const
    {
      return NonzeroColumnRange(NonzeroIterator<false>(*this, column),
        NonzeroIterator<false>(*this, column, 0));
    }

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

    void pivot2(Index row, Index column)
    {
      assert(get(row, column) != 0);

      for (Index i = 0; i < numRows(); ++i)
      {
        if (i == row)
          continue;

        const Value& first = get(i, column);
        if (first == 0)
          continue;

        for (Index j = 0; j < numColumns(); ++j)
        {
          if (j == column)
            continue;

          const Value& second = get(row, j);
          if (second == 0)
            continue;

          set(i, j, 1 - get(i, j));
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

    void pivot3(Index row, Index column)
    {
      const Value& base = get(row, column);
      assert(base != 0);
      signed char denominator = base > 0 ? 1 : -1;

      for (Index i = 0; i < numRows(); ++i)
      {
        if (i == row)
          continue;

        const Value& first = get(i, column);
        if (first == 0)
          continue;
        signed char firstInt = first > 0 ? 1 : -1;

        for (Index j = 0; j < numColumns(); ++j)
        {
          if (j == column)
            continue;

          const Value& second = get(row, j);
          if (second == 0)
            continue;
          signed char numerator = firstInt * (second > 0 ? 1 : -1);

          const Value& value = get(i, j);
          signed char result = (int(get(i,j)) - numerator / denominator + 4) % 3 - 1;
          assert(result >= -1 && result <= 1);
          set(i, j,  result);
        }
      }
    }

  protected:
    Value* _data; /// Matrix entries with row-major ordering.
    Index _numRows; /// Number of rows.
    Index _numColumns; /// Number of columns.
  };


  /**
   * Exception to indicate a pivot on a zero element.
   */

  class matrix_binary_pivot_exception: public std::exception
  {
  public:
    const char* what() const throw ()
    {
      return "Cannot pivot on a zero entry!";
    }
  };

  /**
   * Free function to set a matrix value of a permuted matrix.
   *
   * @param matrix The permuted matrix
   * @param row Row index
   * @param column Column index
   * @param value New value
   */

  template <typename MatrixType>
  inline void matrix_set_value(MatrixType& matrix, size_t row, size_t column, typename MatrixType::value_type value)
  {
    matrix(row, column) = value;
  }

  /**
   * Dummy function for the version with writable orignal matrix.
   */

  template <typename MatrixType>
  inline void matrix_set_value(const MatrixType& matrix, size_t row, size_t column, typename MatrixType::value_type value)
  {
    assert (false);
  }

  /**
   * Free function to permute two rows of a permuted matrix.
   *
   * @param matrix The permuted matrix
   * @param index1 First index
   * @param index2 Second index
   */

  template <typename MatrixType>
  inline void matrix_permute1(MatrixType& matrix, size_t index1, size_t index2)
  {
    for (size_t index = 0; index < matrix.size2(); ++index)
    {
      std::swap(matrix(index1, index), matrix(index2, index));
    }
  }

  /**
   * Free function to permute two columns of a permuted matrix.
   *
   * @param matrix The permuted matrix
   * @param index1 First index
   * @param index2 Second index
   */

  template <typename MatrixType>
  inline void matrix_permute2(MatrixType& matrix, size_t index1, size_t index2)
  {
    for (size_t index = 0; index < matrix.size1(); ++index)
    {
      std::swap(matrix(index, index1), matrix(index, index2));
    }
  }

  /**
   * Free function to perform a binary pivot on a matrix.
   *
   * @param matrix The matrix
   * @param i Row index
   * @param j Column index
   */

  template <typename MatrixType>
  void matrix_binary_pivot(MatrixType& matrix, size_t i, size_t j)
  {
    typedef typename MatrixType::value_type value_type;
    const value_type& base_value = matrix(i, j);

    if (base_value == 0)
    {
      throw matrix_binary_pivot_exception();
    }

    for (size_t row = 0; row < matrix.size1(); ++row)
    {
      if (row == i)
      {
        continue;
      }
      const value_type& first = matrix(row, j);
      if (first == 0)
      {
        continue;
      }

      for (size_t column = 0; column < matrix.size2(); ++column)
      {
        if (column == j)
        {
          continue;
        }
        const value_type& second = matrix(i, column);
        if (second == 0)
        {
          continue;
        }
        matrix(row, column) = 1 - matrix(row, column);
      }
    }
  }

  /**
   * Free function to perform a ternary pivot on a matrix.
   *
   * @param matrix The matrix
   * @param i Row index
   * @param j Column index
   */

  template <typename MatrixType>
  void matrix_ternary_pivot(MatrixType& matrix, size_t i, size_t j)
  {
    typedef typename MatrixType::value_type value_type;
    const value_type& base_value = matrix(i, j);

    if (base_value == 0)
    {
      throw matrix_binary_pivot_exception();
    }

    for (size_t row = 0; row < matrix.size1(); ++row)
    {
      if (row == i)
      {
        continue;
      }
      const value_type& first = matrix(row, j);
      if (first == 0)
      {
        continue;
      }

      for (size_t column = 0; column < matrix.size2(); ++column)
      {
        if (column == j)
        {
          continue;
        }
        const value_type& second = matrix(i, column);
        if (second == 0)
        {
          continue;
        }
        value_type value = matrix(row, column) - first * second / base_value;
        while (value > 1)
          value -= 3;
        while (value < -1)
          value += 3;
        matrix(row, column) = value;
      }
    }
  }

  /**
   * Counts how many of the rows in the specified range have a property.
   *
   * @param matrix The given matrix
   * @param row_first First row in range
   * @param row_beyond Beyond row in range
   * @param column_first First column for property check
   * @param column_beyond Beyond column for property check
   * @param check Property checking routine
   * @return Number of rows with the property
   */

  template <typename MatrixType, typename PropertyCheck>
  inline size_t matrix_count_property_row_series(const MatrixType& matrix, size_t row_first, size_t row_beyond, size_t column_first,
      size_t column_beyond, PropertyCheck check)
  {
    for (size_t row = row_first; row < row_beyond; ++row)
    {
      for (size_t column = column_first; column < column_beyond; ++column)
        check(matrix(row, column));
      if (!check())
        return row - row_first;
    }
    return row_beyond - row_first;
  }

  /**
   * Counts how many of the columns in the specified range have a property.
   *
   * @param matrix The given matrix
   * @param row_first First column for property check
   * @param row_beyond Beyond column for property check
   * @param column_first First row in range
   * @param column_beyond Beyond row in range
   * @param check Property checking routine
   * @return Number of rows with the property
   */

  template <typename MatrixType, typename PropertyCheck>
  inline size_t matrix_count_property_column_series(const MatrixType& matrix, size_t row_first, size_t row_beyond, size_t column_first,
      size_t column_beyond, PropertyCheck check)
  {
    matrix_transposed <const MatrixType> transposed(matrix);
    return matrix_count_property_row_series(transposed, column_first, column_beyond, row_first, row_beyond, check);
  }

  /**
   * Prints a matrix
   *
   * @param matrix The matrix to be printed
   */

  template <typename MatrixType>
  inline void matrix_print(const MatrixType& matrix)
  {
    for (size_t row = 0; row < matrix.size1(); ++row)
    {
      for (size_t column = 0; column < matrix.size2(); ++column)
      {
        std::cout << " " << std::setw(2) << matrix(row, column);
      }
      std::cout << "\n";
    }
    std::cout << std::flush;
  }

  /**
   * Tests two matrices for exact equality.
   *
   * @param matrix1 First matrix
   * @param matrix2 Second matrix
   * @return true if and only if all entries match
   */

  template <typename MatrixType1, typename MatrixType2>
  inline bool matrix_equals(const MatrixType1& matrix1, const MatrixType2& matrix2)
  {
    if (matrix1.size1() != matrix2.size1())
      return false;
    if (matrix1.size2() != matrix2.size2())
      return false;

    for (size_t row = 0; row < matrix1.size1(); ++row)
    {
      for (size_t column = 0; column < matrix1.size2(); ++column)
      {
        if (matrix1(row, column) != matrix2(row, column))
          return false;
      }
    }
    return true;
  }

  template <typename Matrix>
  bool matrix_row_zero(const Matrix& matrix, size_t row, size_t column_first, size_t column_beyond)
  {
    for (size_t c = column_first; c != column_beyond; ++c)
      if (matrix(row, c) != 0)
        return false;
    return true;
  }

  template <typename Matrix>
  bool matrix_column_zero(const Matrix& matrix, size_t column, size_t row_first, size_t row_beyond)
  {
    return matrix_row_zero(make_transposed_matrix(matrix), column, row_first, row_beyond);
  }

  template <typename Matrix>
  bool find_smallest_nonzero_matrix_entry(const Matrix& matrix, size_t row_first, size_t row_beyond, size_t column_first, size_t column_beyond,
      size_t& row, size_t& column)
  {
    bool result = false;
    int current_value = 0;
    for (size_t r = row_first; r != row_beyond; ++r)
    {
      for (size_t c = column_first; c != column_beyond; ++c)
      {
        int value = matrix(r, c);
        if (value == 0)
          continue;

        value = value >= 0 ? value : -value;

        if (!result || value < current_value)
        {
          result = true;
          row = r;
          column = c;
          current_value = value;
        }
      }
    }
    return result;
  }

  template <typename Matrix1, typename Matrix2>
  bool equals(const Matrix1& first, const Matrix2& second)
  {
    if (first.size1() != second.size1())
      return false;
    if (first.size2() != second.size2())
      return false;
    for (size_t r = 0; r < first.size1(); ++r)
    {
      for (size_t c = 0; c < first.size2(); ++c)
      {
        if (first(r, c) != second(r, c))
          return false;
      }
    }
    return true;
  }
} /* namespace tu */
