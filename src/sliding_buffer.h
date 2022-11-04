#ifndef AVP_SLIDINGBUFFER_H_
#define AVP_SLIDINGBUFFER_H_

#include <map>
#include <utility>
#include <vector>

namespace avp {
template <typename Time, typename T>
class SlidingBuffer {
 public:
  SlidingBuffer() = delete;
  SlidingBuffer(const size_t capacity, const Time max_interval)
      : max_interval_(max_interval), capacity_(capacity) {}

  ~SlidingBuffer() {}

  void Put(const Time& timestamp, const T& data) {
    if (buf_.size() >= capacity_) {
      if (timestamp > StartTime()) {
        buf_.erase(buf_.begin());
      } else {
        return;
      }
    }
    buf_.emplace(timestamp, data);
  }

  void GetElements(const Time& from, const Time& to,
                   std::vector<std::pair<Time, T>>* elements) const {
    if (buf_.empty()) return;
    typename std::map<Time, T>::const_iterator from_iter =
        buf_.upper_bound(from);
    typename std::map<Time, T>::const_iterator to_iter = buf_.upper_bound(to);

    if (from_iter != buf_.begin()) from_iter--;
    while (from_iter != to_iter) {
      elements->push_back({from_iter->first, from_iter->second});
      from_iter++;
    }
  }

  bool GetNearestElements(const Time& time, T* element) const {
    if (buf_.empty()) {
      return false;
    }
    typename std::map<Time, T>::const_iterator upper = buf_.upper_bound(time);
    typename std::map<Time, T>::const_iterator lower = buf_.lower_bound(time);

    typename std::map<Time, T>::const_iterator nearest_it = buf_.end();
    if (upper == buf_.end()) {
      if (lower == buf_.end()) lower--;
      nearest_it = lower;
    } else {
      if (lower == upper) lower--;
      if (lower == buf_.end()) {
        nearest_it = upper;
      } else {
        if (std::abs(lower->first - time) > std::abs(upper->first - time)) {
          nearest_it = upper;
        } else {
          nearest_it = lower;
        }
      }
    }

    bool ret = false;
    if (std::abs(nearest_it->first - time) < max_interval_) {
      *element = nearest_it->second;
      ret = true;
    }

    return ret;
  }

  Time StartTime() {
    if (buf_.empty()) {
      return Time(-1);
    }
    return buf_.begin()->first;
  }

  Time EndTime() {
    if (buf_.empty()) {
      return Time(-1);
    }
    return buf_.rbegin()->first;
  }

  void Clear() { buf_.clear(); }

  size_t Capacity() { return capacity_; }

  size_t Size() { return buf_.size(); }

  std::map<Time, T>& GetBuffer() { return buf_; }

 private:
  Time max_interval_;
  size_t capacity_;
  std::map<Time, T> buf_;
};

}  // namespace avp

#endif
