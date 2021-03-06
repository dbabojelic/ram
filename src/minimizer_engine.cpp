// Copyright (c) 2020 Robert Vaser

#include "ram/minimizer_engine.hpp"

#include <cassert>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {

static std::uint64_t First(const std::pair<std::uint64_t, std::uint64_t>& pr) {
  return pr.first;
}

static std::uint64_t Second(const std::pair<std::uint64_t, std::uint64_t>& pr) {
  return pr.second;
}

}  // namespace

namespace ram {

// clang-format off
const std::vector<std::uint64_t> kCoder = {
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255,   0, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255,   0,   1 ,  1,   0, 255, 255,   2,
      3, 255, 255,   2, 255,   1,   0, 255,
    255, 255,   0,   1,   3,   3,   2,   0,
    255,   3, 255, 255, 255, 255, 255, 255,
    255,   0,   1,   1,   0, 255, 255,   2,
      3, 255, 255,   2, 255,   1,   0, 255,
    255, 255,   0,   1,   3,   3,   2,   0,
    255,   3, 255, 255, 255, 255, 255, 255};
// clang-format on

MinimizerEngine::MinimizerEngine(
    std::uint32_t kmer_len, std::uint32_t window_len,
    std::uint32_t chaining_score_treshold,
    std::uint64_t chain_enlongation_stop_criteria,
    std::uint8_t chain_minimizer_cnt_treshold, std::uint32_t best_n,
    std::uint32_t reduce_win_sz, bool hpc, bool robust_winnowing,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool)
    : k_(std::min(std::max(kmer_len, 1U), 32U)),
      w_(window_len),
      occurrence_(-1),
      m_(chaining_score_treshold),
      g_(chain_enlongation_stop_criteria),
      n_(chain_minimizer_cnt_treshold),
      best_n_(best_n),
      reduce_win_sz_(reduce_win_sz),
      robust_winnowing_(robust_winnowing),
      hpc_(hpc),
      minimizers_(1U << std::min(14U, 2 * k_)),
      index_(minimizers_.size()),
      thread_pool_(thread_pool ? thread_pool
                               : std::make_shared<thread_pool::ThreadPool>(1)) {
}

void MinimizerEngine::Minimize(
    std::vector<std::unique_ptr<biosoup::Sequence>>::const_iterator begin,
    std::vector<std::unique_ptr<biosoup::Sequence>>::const_iterator end) {
  for (auto& it : minimizers_) {
    it.clear();
  }
  for (auto& it : index_) {
    it.clear();
  }

  if (begin >= end) {
    return;
  }

  {
    std::uint64_t bin_mask = minimizers_.size() - 1;

    std::vector<std::future<std::vector<uint128_t>>> futures;
    for (auto it = begin; it != end; ++it) {
      futures.emplace_back(thread_pool_->Submit(
          [&](std::vector<std::unique_ptr<biosoup::Sequence>>::const_iterator
                  it)  // NOLINT
          -> std::vector<uint128_t> { return Minimize(*it); },
          it));
    }
    for (auto& it : futures) {
      for (const auto& jt : it.get()) {
        minimizers_[jt.first & bin_mask].emplace_back(jt);
      }
    }
  }

  {
    std::vector<std::future<void>> futures;
    for (std::uint32_t i = 0; i < minimizers_.size(); ++i) {
      if (minimizers_[i].empty()) {
        continue;
      }

      futures.emplace_back(thread_pool_->Submit(
          [&](std::uint32_t bin) -> void {
            RadixSort(minimizers_[bin].begin(), minimizers_[bin].end(), k_ * 2,
                      ::First);

            for (std::uint64_t i = 0, c = 0; i < minimizers_[bin].size(); ++i) {
              if (i > 0 && minimizers_[bin][i - 1].first !=
                               minimizers_[bin][i].first) {  // NOLINT
                index_[bin].emplace(minimizers_[bin][i - 1].first,
                                    std::make_pair(i - c, c));
                c = 0;
              }
              if (i == minimizers_[bin].size() - 1) {
                index_[bin].emplace(minimizers_[bin][i].first,
                                    std::make_pair(i - c, c + 1));
              }
              ++c;
            }
          },
          i));
    }
    for (const auto& it : futures) {
      it.wait();
    }
  }
}

void MinimizerEngine::Filter(double frequency) {
  if (!(0 <= frequency && frequency <= 1)) {
    throw std::invalid_argument(
        "[ram::MinimizerEngine::Filter] error: invalid frequency");
  }

  if (frequency == 0) {
    occurrence_ = -1;
    return;
  }

  std::vector<std::uint32_t> occurrences;
  for (const auto& it : index_) {
    for (const auto& jt : it) {
      occurrences.emplace_back(jt.second.second);
    }
  }

  if (occurrences.empty()) {
    occurrence_ = -1;
    return;
  }

  std::nth_element(occurrences.begin(),
                   occurrences.begin() + (1 - frequency) * occurrences.size(),
                   occurrences.end());
  occurrence_ = occurrences[(1 - frequency) * occurrences.size()] + 1;
}

std::vector<biosoup::Overlap> MinimizerEngine::Map(
    const std::unique_ptr<biosoup::Sequence>& sequence, bool avoid_equal,
    bool avoid_symmetric, bool micromize, double micromize_factor,
    std::uint8_t N) const {
  auto sketch = Minimize(sequence, micromize, micromize_factor, N);
  if (sketch.empty()) {
    return std::vector<biosoup::Overlap>{};
  }

  std::uint64_t bin_mask = minimizers_.size() - 1;
  std::vector<uint128_t> matches;
  for (const auto& it : sketch) {
    std::uint32_t bin = it.first & bin_mask;
    auto match = index_[bin].find(it.first);
    if (match == index_[bin].end() || match->second.second > occurrence_) {
      continue;
    }

    auto jt = minimizers_[bin].begin() + match->second.first;
    auto end = jt + match->second.second;
    for (; jt != end; ++jt) {
      std::uint64_t rhs_id = jt->second >> 32;
      if (avoid_equal && static_cast<std::uint64_t>(sequence->id) == rhs_id) {
        continue;
      }
      if (avoid_symmetric &&
          static_cast<std::uint64_t>(sequence->id) > rhs_id) {  // NOLINT
        continue;
      }

      std::uint64_t strand = (it.second & 1) == (jt->second & 1);
      std::uint64_t lhs_pos = it.second << 32 >> 33;
      std::uint64_t rhs_pos = jt->second << 32 >> 33;

      std::uint64_t diagonal =
          !strand ? rhs_pos + lhs_pos : rhs_pos - lhs_pos + (3ULL << 30);

      matches.emplace_back((((rhs_id << 1) | strand) << 32) | diagonal,
                           (lhs_pos << 32) | rhs_pos);
    }
  }

  return Chain(sequence->id, std::move(matches));
}

std::vector<biosoup::Overlap> MinimizerEngine::Map(
    const std::unique_ptr<biosoup::Sequence>& lhs,
    const std::unique_ptr<biosoup::Sequence>& rhs, bool micromize,
    std::uint8_t N) const {
  auto lhs_sketch = Minimize(lhs, micromize, 0., N);
  if (lhs_sketch.empty()) {
    return std::vector<biosoup::Overlap>{};
  }

  auto rhs_sketch = Minimize(rhs);
  if (rhs_sketch.empty()) {
    return std::vector<biosoup::Overlap>{};
  }

  RadixSort(lhs_sketch.begin(), lhs_sketch.end(), k_ * 2, ::First);
  RadixSort(rhs_sketch.begin(), rhs_sketch.end(), k_ * 2, ::First);

  std::uint64_t rhs_id = rhs->id;

  std::vector<uint128_t> matches;
  for (std::uint32_t i = 0, j = 0; i < lhs_sketch.size(); ++i) {
    while (j < rhs_sketch.size()) {
      if (lhs_sketch[i].first < rhs_sketch[j].first) {
        break;
      } else if (lhs_sketch[i].first == rhs_sketch[j].first) {
        for (std::uint32_t k = j; k < rhs_sketch.size(); ++k) {
          if (lhs_sketch[i].first != rhs_sketch[k].first) {
            break;
          }

          std::uint64_t strand =
              (lhs_sketch[i].second & 1) == (rhs_sketch[k].second & 1);
          std::uint64_t lhs_pos = lhs_sketch[i].second << 32 >> 33;
          std::uint64_t rhs_pos = rhs_sketch[k].second << 32 >> 33;
          std::uint64_t diagonal =
              !strand ? rhs_pos + lhs_pos : rhs_pos - lhs_pos + (3ULL << 30);

          matches.emplace_back((((rhs_id << 1) | strand) << 32) | diagonal,
                               (lhs_pos << 32) | rhs_pos);
        }
        break;
      } else {
        ++j;
      }
    }
  }

  return Chain(lhs->id, std::move(matches));
}

std::vector<biosoup::Overlap> MinimizerEngine::Chain(
    std::uint64_t lhs_id, std::vector<uint128_t>&& matches) const {
  RadixSort(matches.begin(), matches.end(), 64, ::First);
  matches.emplace_back(-1, -1);  // stop dummy

  std::vector<uint128_t> intervals;
  for (std::uint64_t i = 1, j = 0; i < matches.size(); ++i) {  // NOLINT
    if (matches[i].first - matches[j].first > 500) {
      if (i - j >= n_) {
        if (!intervals.empty() && intervals.back().second > j) {  // extend
          intervals.back().second = i;
        } else {  // new
          intervals.emplace_back(j, i);
        }
      }
      ++j;
      while (j < i && matches[i].first - matches[j].first > 500) {
        ++j;
      }
    }
  }

  std::vector<biosoup::Overlap> dst;
  for (const auto& it : intervals) {
    std::uint64_t j = it.first;
    std::uint64_t i = it.second;

    if (i - j < n_) {
      continue;
    }

    RadixSort(matches.begin() + j, matches.begin() + i, 64, ::Second);

    std::uint64_t strand = matches[j].first >> 32 & 1;

    std::vector<std::uint64_t> indices;
    if (strand) {                    // same strand
      indices = LongestSubsequence(  // increasing
          matches.begin() + j, matches.begin() + i, std::less<std::uint64_t>());
    } else {                         // different strand
      indices = LongestSubsequence(  // decreasing
          matches.begin() + j, matches.begin() + i,
          std::greater<std::uint64_t>());
    }

    if (indices.size() < n_) {
      continue;
    }

    indices.emplace_back(matches.size() - 1 - j);  // stop dummy from above
    for (std::uint64_t k = 1, l = 0; k < indices.size(); ++k) {
      if ((matches[j + indices[k]].second >> 32) -
              (matches[j + indices[k - 1]].second >> 32) >
          g_) {
        if (k - l < n_) {
          l = k;
          continue;
        }

        std::uint32_t lhs_matches = 0;
        std::uint32_t lhs_begin = 0;
        std::uint32_t lhs_end = 0;
        std::uint32_t rhs_matches = 0;
        std::uint32_t rhs_begin = 0;
        std::uint32_t rhs_end = 0;

        for (std::uint64_t m = l; m < k; ++m) {
          std::uint32_t lhs_pos = matches[j + indices[m]].second >> 32;
          if (lhs_pos > lhs_end) {
            lhs_matches += lhs_end - lhs_begin;
            lhs_begin = lhs_pos;
          }
          lhs_end = lhs_pos + k_;

          std::uint32_t rhs_pos = matches[j + indices[m]].second << 32 >> 32;
          rhs_pos = strand ? rhs_pos : (1U << 31) - (rhs_pos + k_ - 1);
          if (rhs_pos > rhs_end) {
            rhs_matches += rhs_end - rhs_begin;
            rhs_begin = rhs_pos;
          }
          rhs_end = rhs_pos + k_;
        }
        lhs_matches += lhs_end - lhs_begin;
        rhs_matches += rhs_end - rhs_begin;
        if (std::min(lhs_matches, rhs_matches) < m_) {
          l = k;
          continue;
        }

        // clang-format off
        dst.emplace_back(
            lhs_id,
            matches[j + indices[l]].second >> 32,  // lhs_begin
            k_ + (matches[j + indices[k - 1]].second >> 32),  // lhs_end
            matches[j].first >> 33,  // rhs_id
            strand ?  // rhs_begin
                matches[j + indices[l]].second << 32 >> 32 :
                matches[j + indices[k - 1]].second << 32 >> 32,
            k_ + (strand ?  // rhs_end
                matches[j + indices[k - 1]].second << 32 >> 32 :
                matches[j + indices[l]].second << 32 >> 32),
            std::min(lhs_matches, rhs_matches),  // score
            strand);
        // clang-format on

        l = k;
      }
    }
  }

  if (best_n_ && best_n_ < dst.size()) {
    // take only best_n_ overlaps
    std::sort(
        dst.begin(), dst.end(),
        [](const biosoup::Overlap& lhs, const biosoup::Overlap& rhs) -> bool {
          return lhs.score > rhs.score;
        });
    dst.resize(best_n_);
  }
  return dst;
}

std::vector<MinimizerEngine::uint128_t> MinimizerEngine::Minimize(
    const std::unique_ptr<biosoup::Sequence>& sequence, bool micromize,
    double micromize_factor, std::uint8_t N) const {
  if (sequence->data.size() < k_) {
    return std::vector<uint128_t>{};
  }

  std::uint64_t mask = (1ULL << (k_ * 2)) - 1;

  auto hash = [&](std::uint64_t key) -> std::uint64_t {
    key = ((~key) + (key << 21)) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key << 31)) & mask;
    return key;
  };

  std::deque<uint128_t> window;
  auto window_add = [&](std::uint64_t minimizer,
                        std::uint64_t location) -> void {  // NOLINT
    while (!window.empty() && window.back().first > minimizer) {
      window.pop_back();
    }
    window.emplace_back(minimizer, location);
  };
  auto robust_pop = [&]() -> void {
    while ((int)window.size() > 1 && window.front().first == window[1].first) {
      window.pop_front();
    }
  };
  auto window_update = [&](std::uint32_t position) -> void {
    bool popped = false;
    while (!window.empty() && (window.front().second << 32 >> 33) < position) {
      window.pop_front();
      popped = true;
    }
    if (robust_winnowing_ && popped) {
      robust_pop();
    }
  };

  std::uint64_t shift = (k_ - 1) * 2;
  std::uint64_t minimizer = 0;
  std::uint64_t reverse_minimizer = 0;
  std::uint64_t id = static_cast<std::uint64_t>(sequence->id) << 32;
  std::uint64_t is_stored = 1ULL << 63;

  std::vector<uint128_t> dst;

  for (std::uint32_t i = 0, win_span = 0, kmer_span = 0, base_cnt = 0;
       i < sequence->data.size(); ++i, ++win_span, ++kmer_span) {
    std::uint64_t c = kCoder[sequence->data[i]];
    if (c == 255ULL) {
      throw std::invalid_argument(
          "[ram::MinimizerEngine::Minimize] error: invalid character");
    }

    // skip homopoly
    if (hpc_ && i && kCoder[sequence->data[i - 1]] == c) {
      continue;
    }

    // found new char
    base_cnt++;

    // remove last from kmer
    if (base_cnt > k_) {
      kmer_span--;
      if (hpc_) {
        auto last_c = kCoder[sequence->data[i - kmer_span - 1]];
        while (kCoder[sequence->data[i - kmer_span]] == last_c) kmer_span--;
      }
    }

    minimizer = ((minimizer << 2) | c) & mask;
    reverse_minimizer = (reverse_minimizer >> 2) | ((c ^ 3) << shift);
    if (base_cnt >= k_) {
      if (minimizer < reverse_minimizer) {
        window_add(hash(minimizer), (i - (kmer_span)) << 1 | 0);
      } else if (minimizer > reverse_minimizer) {
        window_add(hash(reverse_minimizer), (i - (kmer_span)) << 1 | 1);
      }
    }
    if (base_cnt >= (k_) + (w_ - 1U)) {
      auto stop = window.end();
      if (!window.empty() && robust_winnowing_) stop = window.begin() + 1;
      for (auto it = window.begin(); it != stop; ++it) {
        if (it->first != window.front().first) {
          break;
        }
        if (it->second & is_stored) {
          continue;
        }
        dst.emplace_back(it->first, id | it->second);
        it->second |= is_stored;
      }
      win_span--;
      if (hpc_) {
        auto last_c = kCoder[sequence->data[i - win_span - 1]];
        while (kCoder[sequence->data[i - win_span]] == last_c) win_span--;
      }
      window_update(i - win_span);
    }
  }

  if (micromize) {
    std::uint32_t take = sequence->data.size() / k_;
    if (micromize_factor > 0.) {
      take = (int)dst.size() * micromize_factor;
    }
    if (take < dst.size()) {
      if (2 * N <= dst.size())
        RadixSort(dst.begin() + N, dst.end() - N, k_ * 2, ::First);
      if (N < take)
        dst.insert(dst.begin() + take - N, dst.end() - N, dst.end());
      dst.resize(take);
    }
  }
  if (reduce_win_sz_)
    return Reduce(dst);
  else
    return dst;
}

template <typename T>
void MinimizerEngine::RadixSort(std::vector<uint128_t>::iterator begin,
                                std::vector<uint128_t>::iterator end,
                                std::uint8_t max_bits,
                                T compare) {  //  unary comparison function

  if (begin >= end) {
    return;
  }

  std::vector<MinimizerEngine::uint128_t> dst(end - begin);
  auto dst_begin = dst.begin();
  auto dst_end = dst.end();

  std::uint64_t buckets[0x100]{};  // 256 b
  std::uint8_t shift = 0;
  for (; shift < max_bits; shift += 8) {
    std::uint64_t counts[0x100]{};
    for (auto it = begin; it != end; ++it) {
      ++counts[compare(*it) >> shift & 0xFF];
    }
    for (std::uint64_t i = 0, j = 0; i < 0x100; j += counts[i++]) {
      buckets[i] = j;
    }
    for (auto it = begin; it != end; ++it) {
      *(dst_begin + buckets[compare(*it) >> shift & 0xFF]++) = *it;
    }
    std::swap(dst_begin, begin);
    std::swap(dst_end, end);
  }

  if (shift / 8 & 1) {  // copy the sorted array for odd cases
    for (; begin != end; ++begin, ++dst_begin) {
      *dst_begin = *begin;
    }
  }
}

template <typename T>
std::vector<std::uint64_t> MinimizerEngine::LongestSubsequence(
    std::vector<uint128_t>::const_iterator begin,
    std::vector<uint128_t>::const_iterator end,
    T compare) {  // binary comparison function

  if (begin >= end) {
    return std::vector<std::uint64_t>{};
  }

  std::vector<std::uint64_t> minimal(end - begin + 1, 0);
  std::vector<std::uint64_t> predecessor(end - begin, 0);

  std::uint64_t longest = 0;
  for (auto it = begin; it != end; ++it) {
    std::uint64_t lo = 1, hi = longest;
    while (lo <= hi) {
      std::uint64_t mid = lo + (hi - lo) / 2;
      if (((begin + minimal[mid])->second >> 32) < (it->second >> 32) &&
          compare((begin + minimal[mid])->second << 32 >> 32,
                  it->second << 32 >> 32)) {
        lo = mid + 1;
      } else {
        hi = mid - 1;
      }
    }

    predecessor[it - begin] = minimal[lo - 1];
    minimal[lo] = it - begin;
    longest = std::max(longest, lo);
  }

  std::vector<std::uint64_t> dst;
  for (std::uint64_t i = 0, j = minimal[longest]; i < longest; ++i) {
    dst.emplace_back(j);
    j = predecessor[j];
  }
  std::reverse(dst.begin(), dst.end());

  return dst;
}
uint64_t MinimizerEngine::GetMinimizerIndexSize() const {
  uint64_t ret = 0;
  for (const auto& it : minimizers_) {
    ret += it.size();
  }
  return ret;
}
std::vector<MinimizerEngine::uint128_t> MinimizerEngine::Reduce(
    const std::vector<uint128_t>& dst) const {
  std::uint32_t win_sz = reduce_win_sz_;

  if (dst.empty()) return {};

  if (win_sz > dst.size()) {
    int mini = 0;
    for (int i = 1; i < (int)dst.size(); i++)
      if (dst[i].first < dst[mini].first) mini = i;
    return std::vector<uint128_t>{dst[mini]};
  }

  std::vector<uint128_t> ret;
  std::vector<bool> stored(dst.size(), false);

  std::deque<std::pair<std::uint64_t, std::uint32_t>> window;

  auto window_add = [&](std::uint64_t minimizer,
                        std::uint32_t location) -> void {
    while (!window.empty() && window.back().first > minimizer) {
      window.pop_back();
    }
    window.emplace_back(minimizer, location);
  };
  auto window_update = [&](std::uint32_t position) -> void {
    while (!window.empty() && window.front().second < position) {
      window.pop_front();
    }
  };

  auto collect = [&]() -> void {
    for (auto it = window.begin(); it != window.end(); it++) {
      if (it->first != window.front().first) break;
      if (stored[it->second]) continue;
      stored[it->second] = true;
      ret.push_back(dst[it->second]);
    }
  };

  for (uint32_t i = 0; i < win_sz; i++) {
    window_add(dst[i].first, i);
  }

  for (uint32_t i = win_sz; i < dst.size(); i++) {
    collect();
    window_update(i - win_sz + 1);
    window_add(dst[i].first, i);
  }
  collect();

  return ret;
}
std::vector<biosoup::Overlap> MinimizerEngine::MapBeginEnd(
    const std::unique_ptr<biosoup::Sequence>& sequence, bool avoid_equal,
    bool avoid_symmetric, std::uint32_t K) const {
  auto sequence_size = sequence->data.size();
  if (sequence_size <= 4 * K)
    return Map(sequence, avoid_equal, avoid_symmetric);

  auto begin_seq = std::unique_ptr<biosoup::Sequence>(
      new biosoup::Sequence(sequence->name, sequence->data.substr(0, K)));
  auto end_seq = std::unique_ptr<biosoup::Sequence>(new biosoup::Sequence(
      sequence->name, sequence->data.substr(sequence_size - K, K)));

  auto begin_overlap = Map(begin_seq, avoid_equal, avoid_symmetric);
  auto end_overlap = Map(end_seq, avoid_equal, avoid_symmetric);
  if (begin_overlap.empty() || end_overlap.empty()) return {};

  std::uint64_t min_diff = std::numeric_limits<std::uint64_t>::max();
  int ansi = -1;
  int ansj = -1;

  int max_index_sum = begin_overlap.size() + end_overlap.size() - 2;
  double penalty = 1.0;
  const double penalty_mult = 1.08;

  for (int index_sum = 0; index_sum <= max_index_sum; index_sum++) {
    //    bool found = false;

    for (std::uint32_t i = 0, j = index_sum - i;
         j >= 0 && i < begin_overlap.size(); i++, j--) {
      if (j >= end_overlap.size()) {
        continue;
      }

      const auto& bov = begin_overlap[i];
      const auto& eov = end_overlap[j];

      if (bov.strand != eov.strand) continue;
      if (bov.rhs_id != eov.rhs_id) continue;
      auto rhs_begin = bov.rhs_begin;
      auto rhs_end = eov.rhs_end;
      if (!eov.strand) {
        rhs_begin = eov.rhs_begin;
        rhs_end = bov.rhs_end;
      }

      if (rhs_begin > rhs_end) continue;
      int candidate_len = rhs_end - rhs_begin;
      //      if (candidate_len > 0.95 * sequence_size &&
      //          candidate_len < 1.05 * sequence_size) {
      //        ansi = i;
      //        ansj = j;
      //        found = true;
      //        break;
      //      }
      int candi_diff =
          penalty * std::abs(candidate_len - static_cast<int>(sequence_size));
      if (candi_diff < min_diff) {
        ansi = i;
        ansj = j;
        min_diff = candi_diff;
      }
    }
    penalty *= penalty_mult;
  }

  if (ansi == -1) return {};

  auto lhs_id = sequence->id;
  auto lhs_begin = begin_overlap[ansi].lhs_begin;
  std::uint32_t lhs_end = end_overlap[ansj].lhs_end + sequence_size - K;
  auto rhs_id = begin_overlap[ansi].rhs_id;
  auto rhs_begin = begin_overlap[ansi].rhs_begin;
  auto rhs_end = end_overlap[ansj].rhs_end;

  if (!begin_overlap[ansi].strand) {
    lhs_begin = end_overlap[ansj].lhs_begin;
    lhs_end = begin_overlap[ansi].lhs_end + sequence_size - K;
    rhs_begin = end_overlap[ansj].rhs_begin;
    rhs_end = begin_overlap[ansi].rhs_end;
  }

  return {biosoup::Overlap(lhs_id, lhs_begin, lhs_end, rhs_id, rhs_begin,
                           rhs_end,
                           std::max(lhs_end - lhs_begin, rhs_end - rhs_begin),
                           begin_overlap[ansi].strand)};
}

}  // namespace ram
