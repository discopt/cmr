#include "separation.hpp"

#include <cassert>

#include "total_unimodularity.hpp"
#include "matroid.hpp"

namespace tu
{

  /**
   * Constructs a separation which is none.
   */

  separation::separation() :
    split_(0, 0), upper_right_rank_(-1), lower_left_rank_(-1), special_swap_(0, 0)
  {

  }

  /**
   * Constructs a 1-separation.
   *
   * @param split Size of the first component, which must be at upper-left
   */

  separation::separation(split_type split) :
    split_(split), upper_right_rank_(0), lower_left_rank_(0), special_swap_(0, 0)
  {

  }

  /**
   * Constructs a 2-separation.
   *
   * @param split Size of the first component, which must be at upper-left
   * @param witness1 A witnessing one in the rank 1 part
   */

  separation::separation(split_type split, witness_type witness1) :
    split_(split), special_swap_(0, 0)
  {
    witnesses_.push_back(witness1);
    if (witness1.first >= split.first)
    {
      assert (witness1.second < split.second);

      lower_left_rank_ = 1;
      upper_right_rank_ = 0;
    }
    else
    {
      assert (witness1.second >= split.second);

      lower_left_rank_ = 0;
      upper_right_rank_ = 1;
    }
  }

  /**
   * Constructs a 3-separation.
   *
   * @param split Size of the first component, which must be at upper-left
   * @param witness1 First witnessing one in the rank 2 part
   * @param witness2 Second witnessing one in the rank 2 part
   */

  separation::separation(split_type split, witness_type witness1, witness_type witness2) :
    split_(split), special_swap_(0, 0)
  {
    witnesses_.push_back(witness1);
    witnesses_.push_back(witness2);

    assert (witness1.first != witness2.first);
    assert (witness1.second != witness2.second);

    if (witness1.first >= split.first)
    {
      assert (witness1.second < split.second);

      if (witness2.first >= split.first)
      {
        assert (witness2.second < split.second);

        lower_left_rank_ = 2;
        upper_right_rank_ = 0;
      }
      else
      {
        assert (witness2.second >= split.second);

        lower_left_rank_ = 1;
        upper_right_rank_ = 1;
      }
    }
    else
    {
      assert (witness1.second >= split.second);

      if (witness2.first >= split.first)
      {
        assert (witness2.second < split.second);

        lower_left_rank_ = 1;
        upper_right_rank_ = 1;
      }
      else
      {
        assert (witness2.second >= split.second);

        lower_left_rank_ = 0;
        upper_right_rank_ = 2;
      }
    }
  }

  /**
   * Copy constructor
   *
   * @param other Another separation
   */

  separation::separation(const separation& other) :
    split_(other.split_), witnesses_(other.witnesses_), upper_right_rank_(other.upper_right_rank_), lower_left_rank_(other.lower_left_rank_),
        special_swap_(other.special_swap_)
  {

  }

  /**
   * Assignment operator
   *
   * @param other Another separation
   * @return This separation
   */

  separation& separation::operator=(const separation& other)
  {
    split_ = other.split_;
    witnesses_ = other.witnesses_;
    upper_right_rank_ = other.upper_right_rank_;
    lower_left_rank_ = other.lower_left_rank_;
    special_swap_ = other.special_swap_;

    return *this;
  }

  /**
   * Destructor
   */

  separation::~separation()
  {

  }

  /**
   * @return A separation with the transposed details
   */

  separation separation::transposed()
  {
    separation result(*this);
    std::swap(result.split_.first, result.split_.second);
    for (size_t i = 0; i < witnesses_.size(); ++i)
      std::swap(result.witnesses_[i].first, result.witnesses_[i].second);
    std::swap(result.lower_left_rank_, result.upper_right_rank_);
    if (has_special_swap())
    {
      if (has_special_row_swap())
        result.special_swap_.first = 'c';
      else
        result.special_swap_.first = 'r';
    }
    return result;
  }

} /* namespace tu */
