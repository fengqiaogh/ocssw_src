#ifndef FOCS_FILTERITERATOR
#define FOCS_FILTERITERATOR

#include <iterator>
#include <string>
#include <unordered_map>
#include <vector>

namespace focs {
    /*!
        \class focs::FilterIterator

        \brief Iterator for filtering the results of other iterators.

        This class shouldn't be used by external users.  Its API is very odd and
        cumbersome and should/will probably be made better.  As of writing this
        documentation, it's only used for `focs::KvStore` for filtering values
        via group.

        \tparam InputIterator The type of iterator to walk through.
        \tparam UnaryPredicate Function with which to filter values.

        \section Example
        \snippet tests/filter_iterator.cpp FilterIterator example
    */
template <class InputIterator, class UnaryPredicate>
class FilterIterator {
public:
    /*! \brief Sole and cumbersome constructor. 

        \param first Collection's `begin()` iterator.
        \param last Collection's `end()` iterator.
        \param pred Function with which to filter, generally a lambda.  It
         should take a single value of the type stored in the collection and
         return a boolean.
     */
    FilterIterator(InputIterator first, InputIterator last, UnaryPredicate pred) : first_{first}, it_{first}, last_{last}, pred_{pred} {
        while (it_ != last_) {
            if (pred_(*it_)) {
                break;
            }
            it_++;
        }
    }

    /*! \brief Check if the underlying iterators are not equal. */
    bool operator!=(const FilterIterator& it) {
        return it.it_ != this->it_;
    }

    /*! \brief Check if the underlying iterators are equal. */
    bool operator==(const FilterIterator& it) {
        return it.it_ == this->it_;
    }

    /*! \brief Return the underlying iterator. */
    InputIterator operator*() {
        return it_;
    }

    /*! \brief Find the next value that fulfills the UnaryPredicate's condition, or return end(). */
    FilterIterator operator++() {
        while (it_ != last_) {
            it_++;
            if (it_ == last_ || pred_(*it_)) {
                return *this;
            }
        }
        return end();
    }

    /*! \brief Return a copy of the iterator reset back to the start. */
    FilterIterator begin() const {return FilterIterator(first_, last_, pred_);}
    /*! \brief Return an iterator with no values left to check. */
    FilterIterator end() const {return FilterIterator(last_, last_, pred_);}
private:
    const InputIterator first_; /*!< The beginning of the input collection. */
    InputIterator it_; /*!< The current position of the input collection. */
    const InputIterator last_; /*!< The end of the input collection. */
    const UnaryPredicate pred_; /*!< The filtering function. */
};
} // namespace focs

#endif // FOCS_FILTERITERATOR
