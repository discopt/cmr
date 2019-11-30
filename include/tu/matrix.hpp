#pragma once

#include <tu/config.h>
#include <tu/export.h>

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
      const Value& value;

      Nonzero(Index r, Index c, const Value& v)
        : row(r), column(c), value(v)
      {

      }
    };
  };

  /**
   * \brief Dense matrix with entries of type \p V.
   */

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
   * \brief Sparse matrix with entries of type \p V.
   */

  template<typename V>
  class SparseMatrix : public Matrix<V>
  {
  public:
    typedef typename Matrix<V>::Value Value;
    typedef typename Matrix<V>::Index Index;
    typedef typename Matrix<V>::Nonzero Nonzero;

    /**
     * \brief Matrix data, row-wise or column-wise.
     *
     * Matrix data, row-wise or column-wise. The term major refers to the rows if it is a row-wise
     * view and to columns otherwise. The term minor refers to the opposite.
     */

    struct Data
    {
      /// Array that maps a major to its nonzero range specified by first and beyond entries.
      std::vector<std::pair<Index, Index>> range;
      /// Array that maps the nonzero indices to pairs of minor index and value.
      std::vector<std::pair<Index, Value>> entries;
      /// Whether the data is sorted by minor.
      bool isSorted;

      /**
       * \brief Default constructor.
       */

      Data()
        : isSorted(true)
      {

      }

      /**
       * \brief Move constructor.
       */

      Data(Data&& other)
        : range(std::move(other.range)), entries(std::move(other.entries)), isSorted(other.isSorted)
      {

      }

      /**
       * \brief Copy constructor.
       */

      Data(const Data& other)
        : range(other.range), entries(other.entries), isSorted(other.isSorted)
      {

      }

      /**
       * \brief Move operator.
       */

      Data& operator=(Data&& other)
      {
        range = std::move(other.range);
        entries = std::move(other.entries);
        isSorted = other.isSorted;
        return *this;
      }

      /**
       * \brief Assignment operator.
       */

      Data& operator=(const Data& other)
      {
        range = other.range;
        entries = other.entries;
        isSorted = other.isSorted;
        return *this;
      }

      /**
       * \brief Destructor.
       */

      ~Data() = default;

      /**
       * \brief Computes the transpose version of \p data. \p numMinor is the new number of minor
       * indices.
       */

      void constructFromTransposed(const Data& transposedData, Index numMinor)
      {
        // Iterate over entries and count how many in each major.
        range.resize(numMinor);
        for (auto& r : range)
          r.second = 0;

        for (Index i = 0; i < Index(transposedData.range.size()); ++i)
        {
          for (Index j = transposedData.range[i].first; j < transposedData.range[i].second; ++j)
          {
            ++range[transposedData.entries[j].first].second;
          }
        }

        // Initialize start of each major. This will be incremented when copying.
        range[0].first = 0;
        for (Index i = 1; i < Index(range.size()); ++i)
        {
          range[i].first = range[i-1].first + range[i-1].second;
          range[i-1].second = range[i].first;
        }
        range.back().second = transposedData.entries.size();

        // Copy entries major-wise into slots indicated by first[i], incrementing the latter.
        entries.resize(transposedData.entries.size());
        for (Index transposedMajor = 0; transposedMajor < Index(transposedData.range.size());
          ++transposedMajor)
        {
          for (Index j = transposedData.range[transposedMajor].first;
            j < transposedData.range[transposedMajor].second; ++j)
          {
            auto transposedEntry = transposedData.entries[j];
            entries[range[transposedEntry.first].first] = std::make_pair(transposedMajor,
              transposedEntry.second);
            ++range[transposedEntry.first].first;
          }
        }

        // Re-initialize start of each major.
        range[0].first = 0;
        for (Index i = 1; i < Index(range.size()); ++i)
          range[i].first = range[i-1].second;

        isSorted = true;
      }

      /**
       * \brief Checks representation for consistency, raising a \c std::runtime_error if
       * inconsistent.
       */

      void ensureConsistency(Index numMinor = std::numeric_limits<Index>::max()) const
      {
        for (Index i = 0; i < Index(range.size()); ++i)
        {
          if (range[i].first > entries.size())
            throw std::runtime_error(
                "Inconsistent SparseMatrix::Data: range.first contains index out of range.");
          if (range[i].second > entries.size())
            throw std::runtime_error(
                "Inconsistent SparseMatrix::Data: range.second contains index out of range.");
          if (range[i].first > range[i].second)
            throw std::runtime_error(
                "Inconsistent SparseMatrix::Data: range.first > range.beyond.");
        }
        for (Index i = 0; i < Index(entries.size()); ++i)
        {
          if (entries[i].first >= numMinor)
            throw std::runtime_error("Inconsistent SparseMatrix::Data: minor entry exceeds bound.");
          if (entries[i].second == 0)
            throw std::runtime_error("Inconsistent SparseMatrix::Data: zero entry found.");
        }
        for (Index i = 0; i < Index(range.size()); ++i)
        {
          for (Index j = range[i].first + 1; j < range[i].second; ++j)
          {
            if (entries[j-1].first == entries[j].first)
            {
              std::stringstream ss;
              ss << "Inconsistent SparseMatrix::Data: minor entries of major " << i
                << " contain duplicate indices " << (j-1) << "->" << entries[j-1].first << " and "
                << j << "->" << entries[j].first << ".";
              throw std::runtime_error(ss.str());
            }
            else if (entries[j-1].first > entries[j].first && isSorted)
            {
              std::stringstream ss;
              ss << "Inconsistent SparseMatrix::Data: minor entries of major " << i
                << " are not sorted for indices " << (j-1) << "->" << entries[j-1].first << " and "
                << j << "->" << entries[j].first << ".";
              throw std::runtime_error(ss.str());
            }
          }
        }
      }

      /**
       * \brief Swaps data with \p other.
       */

      void swap(Data& other)
      {
        range.swap(other.range);
        entries.swap(other.entries);
        std::swap(isSorted, other.isSorted);
      }

      /**
       * \brief Finds entry \p major,\p minor.
       * 
       * Find entry \p major,\p minor. If it exists, the index of the entry is returned, and
       * \c std::numeric_limits<std::size_t>::max() otherwise.
       */

      std::size_t find(Index major, Index minor) const
      {
        if (isSorted)
        {
          // Binary search on interval [lower, upper)

          std::size_t lower = range[major].first;
          std::size_t upper = range[major].second;
          std::size_t mid;
          while (lower < upper)
          {
            mid = (lower + upper) / 2;
            if (entries[i].first < minor)
              lower = mid;
            else if (entries[i].first > minor)
              upper = mid;
            else
              return i;
          }
        }
        else
        {
          // Linear scan.

          for (std::size_t i = range[major].first; i < range[major].second; ++i)
          {
            if (entries[i].first == minor)
              return i;
          }
        }

        return std::numeric_limits<std::size_t>::max();
      }
    };

    /**
     * \brief Iterator for row- or column-wise iteration over the entries.
     */

    template <bool Row>
    struct NonzeroIterator
    {
    public:
      NonzeroIterator(const Data& data, Index major)
        : _data(data), _major(major), _index(data.range[major].first)
      {

      }

      NonzeroIterator(const Data& data, Index major, int dummy)
        : _data(data), _major(major), _index(data.range[major].second)
      {

      }

      /**
       * \brief Copy constructor.
       */

      NonzeroIterator(const NonzeroIterator& other)
        : _data(other._data), _major(other._major), _index(other._index)
      {

      }

      /**
       * \brief Returns the current entry.
       */

      inline const Nonzero operator*() const
      {
        // Statically known, so one branch will be optimized away.
        if (Row)
          return Nonzero(_major, _data.entries[_index].first, _data.entries[_index].second);
        else
          return Nonzero(_data.entries[_index].first, _major, _data.entries[_index].second);
      }

      /**
       * \brief Advances to the next entry in the same major.
       */

      inline void operator++()
      {
        ++_index;
        assert(_index <= _data.range[_major].second);
      }

      /**
       * \brief Compares two iterators for being not equal.
       */

      inline bool operator!=(const NonzeroIterator<Row>& other)
      {
        assert(&_data == &other._data);
        return _major != other._major || _index != other._index;
      }

      private:
        /// Reference to the matrix data.
        const Data& _data;
        /// Major index.
        Index _major;
        /// Index of the entry relative to this major.
        Value _index;
    };

    typedef Range<NonzeroIterator<true>> NonzeroRowRange;
    typedef Range<NonzeroIterator<false>> NonzeroColumnRange;

    /**
     * \brief Constructs a 0x0 matrix.
     */

    SparseMatrix()
      : _zero(0)
    {

    }

    /**
     * \brief Move constructor.
     */

    SparseMatrix(SparseMatrix&& other)
      : _rowData(std::move(other._rowData)), _columnData(std::move(other._columnData)), _zero(0)
    {

    }

    /**
     * \brief Move assignment operator.
     */

    SparseMatrix& operator=(SparseMatrix&& other)
    {
      _rowData = std::move(other._rowData);
      _columnData = std::move(other._columnData);
      return *this;
    }

    /**
     * \brief Constructs a matrix.
     *
     * Constructs a matrix. Given entries may contain zeros, but they are filtered.
     *
     * \param majorRow Indicates if data is row-wise.
     * \param numMajor Number of rows (resp. columns) of the matrix.
     * \param numMinor Number of columns (resp. rows) of the matrix.
     * \param numEntries Number of provided entries.
     * \param first Array of size \p numMajor with first index of each row (resp. column).
     * \param beyond Array of size \p numMajor with beyond index of each row (resp. column).
     *   May be \c NULL in which case the first-entry of the next major is considered, i.e.,
     *   the rows (resp. columns) need to be given in an ordered way.
     * \param entryMinors Array of size \p numEntries with column (resp. row) of each entry.
     * \param entryValues Array of size \p numEntries with value of each entry.
     * \param mayContainZeros Whether we have to check for zero entries.
     * \param isSorted Whether the entries of each row (resp. column) are already sorted.
     */

    SparseMatrix(bool majorRow, Index numMajor, Index numMinor, Index numEntries, Index* first,
      Index* beyond, Index* entryMinors, Value* entryValues, bool mayContainZeros = true,
      bool isSorted = false)
      : _zero(0)
    {
      set(majorRow, numMajor, numMinor, numEntries, first, beyond, entryMinors, entryValues,
        mayContainZeros, isSorted);
    }

    /**
     * \brief Constructs a matrix.
     *
     * Constructs a matrix from a prepared data structure.
     *
     * \param majorRow Indicates if data is row-wise.
     * \param majorIsSorted Whether the data of each major is sorted by minor.
     * \param data Matrix data.
     */

    SparseMatrix(bool majorRow, const Data& data, Index numMinor, bool majorIsSorted = false)
      : _zero(0)
    {
      set(majorRow, data, numMinor, majorIsSorted);
    }

    /**
     * \brief Constructs a matrix.
     *
     * Constructs a matrix from a prepared data structure.
     *
     * \param majorRow Indicates if data is row-wise.
     * \param data Matrix data.
     */

    SparseMatrix(bool majorRow, Data&& data, Index numMinor)
      : _zero(0)
    {
      set(majorRow, data, numMinor);
    }

    /**
     * \brief Constructs a matrix.
     *
     * Constructs a matrix from two prepared data structures.
     *
     * \param rowData Matrix row data.
     * \param columnData Matrix column data.
     */

    SparseMatrix(const Data& rowData, const Data& columnData)
      : _zero(0)
    {
      set(rowData, columnData);
    }

    /**
     * \brief Constructs a matrix.
     *
     * Constructs a matrix from a two prepared data structures.
     *
     * \param rowData Matrix row data.
     * \param columnData Matrix column data.
     */

    SparseMatrix(Data&& rowData, Data&& columnData)
      : _zero(0)
    {
      set(std::move(rowData), std::move(columnData));
    }

    /**
     * \brief Sets the contents of the matrix.
     *
     * Sets the contents of the matrix. Given entries may contain zeros, but they are filtered.
     *
     * \param majorRow Indicates if data is row-wise.
     * \param numMajor Number of rows (resp. columns) of the matrix.
     * \param numMinor Number of columns (resp. rows) of the matrix.
     * \param numEntries Number of provided entries.
     * \param first Array of size \p numMajor with first index of each row (resp. column).
     * \param beyond Array of size \p numMajor with beyond index of each row (resp. column).
     *   May be \c nullptr in which case the first-entry of the next major is considered, i.e.,
     *   the rows (resp. columns) need to be ordered.
     * \param entryMinors Array of size \p numEntries with column (resp. row) of each entry.
     * \param entryValues Array of size \p numEntries with value of each entry.
     * \param mayContainZeros Whether zero entries shall be skipped.
     * \param isSorted Whether the entries of each major are sorted by minor.
     */

    void set(bool majorRow, Index numMajor, Index numMinor, Index numEntries, const Index* first,
      const Index* beyond, const Index* entryMinors, const Value* entryValues,
      bool mayContainZeros = true, bool isSorted = false)
    {
      if (mayContainZeros)
      {
        // Use first run over entries to count nonzeros in each row (resp. column).
        _rowData.range.resize(numMajor);
        for (auto& range : _rowData.range)
          range.first = 0;
        Index current = 0;
        for (Index i = 0; i < numMajor; ++i)
        {
          _rowData.range[i].first = current;
          Index b;
          if (beyond != nullptr)
            b = beyond[i];
          else if (i + 1 < numMajor)
            b = first[i+1];
          else
            b = numEntries;
          for (Index j = first[i]; j < b; ++j)
          {
            if (entryValues[j] != 0)
              ++current;
          }
          _rowData.range[i].second = current;
        }

        // Use second run to copy the nonzeros.
        _rowData.entries.resize(current);
        current = 0;
        for (Index i = 0; i < numMajor; ++i)
        {
          Index b;
          if (beyond != nullptr)
            b = beyond[i];
          else if (i + 1 < numMajor)
            b = first[i+1];
          else
            b = numEntries;
          for (Index j = first[i]; j < b; ++j)
          {
            if (entryValues[j] != 0)
            {
              _rowData.entries[current] = std::make_pair(entryMinors[j], entryValues[j]);
              ++current;
            }
          }
        }
      }
      else
      {
        // If we can trust the input then we can copy directly.
        _rowData.range.resize(numMajor);
        if (beyond != nullptr)
        {
          for (Index i = 0; i < numMajor; ++i)
            _rowData.range[i] = std::make_pair(first[i], beyond[i]);
        }
        else
        {
          for (Index i = 0; i + 1 < numMajor; ++i)
            _rowData.range[i] = std::make_pair(first[i], first[i+1]);
          _rowData.range.back() = std::make_pair(first[numMajor - 1], Index(numEntries));
        }
        _rowData.entries.resize(numEntries);
        for (Index i = 0; i < numEntries; ++i)
          _rowData.entries[i] = std::make_pair(entryMinors[i], entryValues[i]);
      }

      _rowData.isSorted = isSorted;

// TODO: remove
//       if (!isSorted)
//       {
//         // We have to sort each row in the row-wise representation.
//         for (Index row = 0; row < Index(_rowData.range.size()); ++row)
//         {
//           auto compare = [] (const std::pair<Index, Value>& a,
//             const std::pair<Index, Value>& b)
//             {
//               return a.first < b.first;
//             };
// 
//           std::sort(_rowData.entries.begin() + _rowData.range[row].first,
//             _rowData.entries.begin() + _rowData.range[row].second, compare);
//         }
//       }

      // Construct column-wise (resp. row-wise) representation.
      _columnData.constructFromTransposed(_rowData, numMinor);

      if (!majorRow)
        _rowData.swap(_columnData);

  #if !defined(NDEBUG)
      ensureConsistency();
  #endif /* !NDEBUG */
    }

    /**
     * \brief Sets the contents of the matrix.
     *
     * Sets the contents of the matrix from a prepared data structure.
     *
     * \param majorRow Indicates if data is row-wise.
     * \param data Matrix data.
     * \param numMinor Number of columns (resp. rows).
     */

    void set(bool majorRow, const Data& data, Index numMinor)
    {
      // We can trust the input and thus we can copy directly.
      _rowData.range = data.range;
      _rowData.entries = data.entries;
      _rowData.isSorted = data.isSorted;

      // Construct column-wise (resp. row-wise) representation.
      _columnData.constructFromTransposed(_rowData, numMinor);

      if (!majorRow)
        _rowData.swap(_columnData);

  #if !defined(NDEBUG)
      ensureConsistency();
  #endif /* !NDEBUG */
    }

    /**
     * \brief Sets the contents of the matrix.
     *
     * Sets the contents of the matrix from a prepared data structure.
     *
     * \param majorRow Indicates if data is row-wise.
     * \param data Matrix data.
     * \param numMinor Number of columns (resp. rows).
     */

    void set(bool majorRow, Data&& data, Index numMinor)
    {
      _rowData = std::move(data);
      _columnData.constructFromTransposed(_rowData, numMinor);
      if (!majorRow)
        _rowData.swap(_columnData);

#if !defined(NDEBUG)
      ensureConsistency();
#endif /* !NDEBUG */
    }

    /**
     * \brief Sets the contents of the matrix.
     *
     * Sets the contents of the matrix from a prepared data structure.
     *
     * \param rowData Matrix row data.
     * \param columnData Matrix column data.
     */

    void set(const Data& rowData, const Data& columnData)
    {
      _rowData.range = rowData.range;
      _rowData.entries = rowData.entries;
      _rowData.isSorted = rowData.isSorted;
      _rowData.ensureConsistency();

      _columnData.range = columnData.range;
      _columnData.entries = columnData.entries;
      _columnData.isSorted = columnData.isSorted;
      _columnData.ensureConsistency();

#if !defined(NDEBUG)
      ensureConsistency();
#endif /* !NDEBUG */
    }

    /**
     * \brief Sets the contents of the matrix.
     *
     * Sets the contents of the matrix from a prepared data structure.
     *
     * \param rowData Matrix row data.
     * \param columnData Matrix column data.
     */

    void set(Data&& rowData, Data&& columnData)
    {
      _rowData = std::move(rowData);
      _columnData = std::move(columnData);
#if !defined(NDEBUG)
      ensureConsistency();
#endif /* !NDEBUG */
    }

    /**
     * \brief Destructor.
     */

    ~SparseMatrix() = default;

    /**
     * \brief Returns the number of rows.
     */

    std::size_t numRows() const
    {
      return _rowData.range.size();
    }

    /**
     * \brief Returns the number of columns.
     */

    std::size_t numColumns() const
    {
      return _columnData.range.size();
    }

    /**
     * \brief Indicates if the row data is sorted by column.
     */

    bool hasSortedRows() const
    {
      return _rowData.isSorted;
    }

    /**
     * \brief Indicates if the column data is sorted by row.
     */

    bool hasSortedColumns() const
    {
      return _columnData.isSorted;
    }

    /**
     * \brief Returns entry at \p row, \p column.
     * 
     * Returns entry at \p row, \p column. If the row or column data is sorted (resp. not sorted)
     * then the time is logarithmic (resp. linear) in the number of nonzeros of the row or column.
     */

    const Value& get(Index row, Index column) const
    {
      assert(row < numRows());
      assert(column < numColumns());

      if (_rowData.isSorted)
      {
        std::size_t columnIndex = _rowData.find(row, column);
        if (columnIndex < std::numeric_limits<std::size_t>::max())
          return _rowData.entries[columnIndex].value;
      }
      else
      {
        std::size_t rowIndex = _columnData.find(column, row);
        if (rowIndex < std::numeric_limits<std::size_t>::max())
          return _columnData.entries[rowIndex].value;
      }

      return _zero;
    }

    /**
     * \brief Returns the number of nonzeros.
     */

    std::size_t numNonzeros() const
    {
      return _rowData.entries.size();
    }

    /**
     * \brief Returns the number of nonzeros in \p row.
     */

    std::size_t countRowNonzeros(Index row) const
    {
      return _rowData.range[row].second - _rowData.range[row].first;
    }

    /**
     * \brief Returns the number of nonzeros in \p column.
     */

    std::size_t countColumnNonzeros(Index column) const
    {
      return _columnData.range[column].second - _columnData.range[column].first;
    }

    /**
     * \brief Returns a range for iterating over the nonzeros of \p row.
     */

    NonzeroRowRange iterateRowNonzeros(Index row) const
    {
      return NonzeroRowRange(NonzeroIterator<true>(_rowData, row),
        NonzeroIterator<true>(_rowData, row, 0));
    }

    /**
     * \brief Returns a range for iterating over the nonzeros of \p column.
     */

    NonzeroColumnRange iterateColumnNonzeros(Index column) const
    {
      return NonzeroColumnRange(NonzeroIterator<false>(_columnData, column),
        NonzeroIterator<false>(_columnData, column, 0));
    }

    /**
     * \brief Transposes the matrix.
     *
     * Transposes the matrix in constant time.
     */

    void transpose()
    {
      _rowData.swap(_columnData);
    }

    /**
     * \brief Transposes the matrix.
     *
     * Transposes the matrix in constant time.
     */

    SparseMatrix transposed() const
    {
      return SparseMatrix(_columnData, _rowData);
    }

    /**
     * \brief Consistency check.
     *
     * Checks consistency of row- and column-data, which includes that all entries are nonzeros
     * and that entries are sorted within each row/column.
     */

    void ensureConsistency() const
    {
      // First ensure individual consistency of row and column data.
      _rowData.ensureConsistency();
      _columnData.ensureConsistency();

      // Copy of a nonzero.
      struct NZ
      {
        Index row;
        Index column;
        Value value;
      };

      // Comparison for entries.
      auto compare = [](const NZ& a, const NZ& b)
      {
        if (a.row != b.row)
          return a.row < b.row;
        if (a.column != b.column)
          return a.column < b.column;
        return a.value < b.value;
      };

      // Construct sorted vector of entries from row data.
      std::vector<NZ> rowNonzeros;
      for (Index i = 0; i < Index(_rowData.range.size()); ++i)
      {
        for (Index j = _rowData.range[i].first; j < _rowData.range[i].second; ++j)
        {
          NZ nz = { i, _rowData.entries[j].first, _rowData.entries[j].second };
          rowNonzeros.push_back(nz);
        }
      }
      std::sort(rowNonzeros.begin(), rowNonzeros.end(), compare);

      // Construct sorted vector of entries from column data.
      std::vector<NZ> columnNonzeros;
      for (Index i = 0; i < Index(_columnData.range.size()); ++i)
      {
        for (Index j = _columnData.range[i].first; j < _columnData.range[i].second; ++j)
        {
          NZ nz = { _columnData.entries[j].first, i, _columnData.entries[j].second };
          columnNonzeros.push_back(nz);
        }
      }
      std::sort(columnNonzeros.begin(), columnNonzeros.end(), compare);

      for (std::size_t i = 0; i < rowNonzeros.size(); ++i)
      {
        if (rowNonzeros[i].row != columnNonzeros[i].row
          || rowNonzeros[i].column != columnNonzeros[i].column
          || rowNonzeros[i].value != columnNonzeros[i].value)
        {
          throw std::runtime_error("Inconsistent Matrix: row and column data differ.");
        }
      }
    }

  private:
    /// Data for row-wise access.
    Data _rowData;
    /// Data for column-wise access.
    Data _columnData;
    /// Zero value.
    Value _zero;
    /// Whether the row and column data is sorted (by columns and rows, respectively).
    bool _isSorted;
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
