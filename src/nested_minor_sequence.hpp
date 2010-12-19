/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef NESTED_MINOR_SEQUENCE_HPP_
#define NESTED_MINOR_SEQUENCE_HPP_

#include <vector>

namespace tu
{

  /**
   * Models a sequence of nested matroid-minors by storing each nesting-step
   * with type and size information.
   */

  class nested_minor_sequence
  {
  public:

    /**
     * Different types of extensions
     */

    enum extension_type
    {
      FIRST_EXTENSION_TYPE = -2,
      ONE_ROW = -2,
      ONE_ROW_TWO_COLUMNS = -1,
      ONE_ROW_ONE_COLUMN = 0,
      TWO_ROWS_ONE_COLUMN = 1,
      ONE_COLUMN = 2,
      BEYOND_EXTENSION_TYPE = 3
    };

    /**
     * Constructs a sequence which only consists of a W3-minor.
     */

    nested_minor_sequence();

    /**
     * Destructor
     */

    virtual ~nested_minor_sequence();

    /**
     * Augments the sequence with a minor of the specified extension type.
     *
     * @param type The extension type of the biggest minor
     */

    void push(extension_type type);

    /**
     * Returns the number of rows in an extension of a given type.
     *
     * @param type Extension type
     * @return Number of rows
     */

    static size_t get_extension_height(extension_type type);

    /**
     * Returns the number of columns in an extension of a given type.
     *
     * @param type Extension type
     * @return Number of columns
     */

    static size_t get_extension_width(extension_type type);

    /**
     * Returns the extension type of a specific minor.
     *
     * @param index Index of the extension. 0 means the extension from W3 to the next
     * @return Extension type
     */

    extension_type get_extension(size_t index) const
    {
      return extensions_[index];
    }

    /**
     * Returns the number of rows in a specific extension
     *
     * @param index Index of the extension. 0 means the extension from W3 to the next
     * @return Number of rows
     */

    size_t get_extension_height(size_t index) const
    {
      return get_extension_height(get_extension(index));
    }

    /**
     * Returns the number of columns in a specific extension
     *
     * @param index Index of the extension. 0 means the extension from W3 to the next
     * @return Number of columns
     */

    size_t get_extension_width(size_t index) const
    {
      return get_extension_width(get_extension(index));
    }

    /**
     * @return Number of rows in the biggest minor
     */

    size_t height() const
    {
      return height_;
    }

    /**
     * @return Number of columns in the biggest minor
     */

    size_t width() const
    {
      return width_;
    }

    /**
     * @return Number of extensions
     */

    size_t size() const
    {
      return extensions_.size();
    }

  private:

    std::vector <extension_type> extensions_;
    size_t height_;
    size_t width_;
  };

  /**
   * Transpose-proxy of a nested minor sequence
   */

  template <typename NestedMinorSequenceType>
  class nested_minor_sequence_transposed
  {
  public:
    /**
     * Constructs the proxy.
     *
     * @param sequence Orinal sequence
     */

    nested_minor_sequence_transposed(NestedMinorSequenceType& sequence) :
      sequence_(sequence)
    {

    }

    /**
     * Copy constructor
     *
     * @param other Another transpose-proxy of a nested minor sequence.
     */

    nested_minor_sequence_transposed(const nested_minor_sequence_transposed <NestedMinorSequenceType>& other) :
      sequence_(other.sequence_)
    {

    }

    /**
     * Destructor
     */

    virtual ~nested_minor_sequence_transposed()
    {

    }

    /**
     * Augments the sequence with a minor of the specified extension type.
     *
     * @param type The extension type of the biggest minor
     */

    void push(nested_minor_sequence::extension_type type)
    {
      sequence_.push(nested_minor_sequence::extension_type(-int(type)));
    }

    /**
     * Returns the extension type of a specific minor.
     *
     * @param index Index of the extension. 0 means the extension from W3 to the next
     * @return Extension type
     */

    nested_minor_sequence::extension_type get_extension(size_t index) const
    {
      return nested_minor_sequence::extension_type(-int(sequence_.get_extension(index)));
    }

    /**
     * Returns the number of rows in a specific extension
     *
     * @param index Index of the extension. 0 means the extension from W3 to the next
     * @return Number of rows
     */

    size_t get_extension_height(size_t index) const
    {
      return sequence_.get_extension_width(index);
    }

    /**
     * Returns the number of columns in a specific extension
     *
     * @param index Index of the extension. 0 means the extension from W3 to the next
     * @return Number of columns
     */

    size_t get_extension_width(size_t index) const
    {
      return sequence_.get_extension_height(index);
    }

    /**
     * @return Number of rows in the biggest minor
     */

    size_t height() const
    {
      return sequence_.width();
    }

    /**
     * @return Number of columns in the biggest minor
     */

    size_t width() const
    {
      return sequence_.height();
    }

    /**
     * @return Number of extensions
     */

    size_t size() const
    {
      return sequence_.size();
    }

  private:
    NestedMinorSequenceType& sequence_;
  };

  /**
   * Creates a transpose proxy of a given nested minor sequence.
   *
   * @param sequence Given sequence
   * @return The transpose proxy
   */

  template <typename NestedMinorSequenceType>
  inline nested_minor_sequence_transposed <NestedMinorSequenceType> make_transposed_nested_minor_sequence(NestedMinorSequenceType& sequence)
  {
    return nested_minor_sequence_transposed <NestedMinorSequenceType> (sequence);
  }

}

#endif /* NESTED_MINOR_SEQUENCE_HPP_ */
