#ifndef SEQ_ITERATOR_HPP
#define SEQ_ITERATOR_HPP

#include "string_view.hpp"
#include <iterator>

class SeqIterator : public std::iterator<std::input_iterator_tag, std::pair<stx::string_view, int>, int>{
  stx::string_view s_;
  bool invalid_ ;
  int lastinvalid_ ;
  int k_ ;

public:
  typedef value_type& reference ;
  typedef value_type* pounter ;
  typedef std::input_iterator_tag iterator_category ;

  typedef int64_t difference_type ;

  SeqIterator(const std::string& s, int k)
    :s_(s), invalid_(true), lastinvalid_(-1), k_(k) {}

  SeqIterator(const SeqIterator& o)
    :s_(o.s), invalid_(o.invalid_), lastinvalid_(o.lastinvalid_), k_(o.k_) {}


