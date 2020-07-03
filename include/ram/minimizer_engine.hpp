// Copyright (c) 2020 Robet Vaser

#ifndef RAM_MINIMIZER_ENGINE_HPP_
#define RAM_MINIMIZER_ENGINE_HPP_

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include "biosoup/overlap.hpp"
#include "biosoup/sequence.hpp"
#include "ram/uint128_t.hpp"
#include "thread_pool/thread_pool.hpp"

namespace ram {

class MinimizerEngine {
 public:
  MinimizerEngine(
      std::uint32_t kmer_len,  // element of [1, 32]
      std::uint32_t window_len,
      std::uint32_t chaining_score_treshold = 100,            // -m param
      std::uint64_t chain_enlongation_stop_criteria = 10000,  // -g param
      std::uint8_t chain_minimizer_cnt_treshold = 4,          // -n param
      std::uint32_t best_n = 0,                               // -b param
      std::uint32_t reduce_win_sz = 0,
      bool robust_winnowing = false,  // -r param
      bool hpc = false,               // use homopolymer-compressed minimizers
      std::shared_ptr<thread_pool::ThreadPool> thread_pool = nullptr);

  MinimizerEngine(const MinimizerEngine&) = delete;
  MinimizerEngine& operator=(const MinimizerEngine&) = delete;

  MinimizerEngine(MinimizerEngine&&) = default;
  MinimizerEngine& operator=(MinimizerEngine&&) = default;

  ~MinimizerEngine() = default;

  // transform set of sequences to minimizer index
  void Minimize(
      std::vector<std::unique_ptr<biosoup::Sequence>>::const_iterator begin,
      std::vector<std::unique_ptr<biosoup::Sequence>>::const_iterator end);

  // set occurrence frequency threshold
  void Filter(double frequency);

  // find overlaps in preconstructed minimizer index
  // micromizers = smallest sequence->data.size() / k minimizers
  std::vector<biosoup::Overlap> Map(
      const std::unique_ptr<biosoup::Sequence>& sequence,
      bool avoid_equal,      // ignore overlaps in which lhs_id == rhs_id
      bool avoid_symmetric,  // ignore overlaps in which lhs_id > rhs_id
      bool micromize = false, std::uint8_t N = 0) const;  // only lhs

  // find overlaps in preconstructed minimizer index
  // using being-end strategy
  std::vector<biosoup::Overlap> MapBeginEnd(
      const std::unique_ptr<biosoup::Sequence>& sequence,
      bool avoid_equal,      // ignore overlaps in which lhs_id == rhs_id
      bool avoid_symmetric,  // ignore overlaps in which lhs_id > rhs_id
      std::uint32_t K) const;

  // find overlaps between a pair of sequences
  std::vector<biosoup::Overlap> Map(
      const std::unique_ptr<biosoup::Sequence>& lhs,
      const std::unique_ptr<biosoup::Sequence>& rhs, bool micromize = false,
      std::uint8_t N = 0) const;  // only lhs

  std::uint64_t GetMinimizerIndexSize() const;

 private:
  using uint64_t_pair = std::pair<std::uint64_t, std::uint64_t>;
  using minimizer_t = std::pair<uint128_t, std::uint64_t>;

  // Match = [127:97] rhs_id
  //         [96:96] strand
  //         [95:64] rhs_pos +- lhs_pos
  //         [63:32] lhs_pos
  //         [31:0] rhs_pos
  std::vector<biosoup::Overlap> Chain(
      std::uint64_t lhs_id, std::vector<uint64_t_pair>&& matches) const;

  // Minimizer = [127:64] kmer
  //             [63:32] id
  //             [31:1] pos
  //             [1:1] strand
  std::vector<minimizer_t> Minimize(
      const std::unique_ptr<biosoup::Sequence>& sequence,
      bool micromize = false, std::uint8_t N = 0) const;

  std::vector<minimizer_t> Reduce(const std::vector<minimizer_t>& dst) const;

  template <typename T>
  static void RadixSort(  // any uint128_t
      std::vector<uint64_t_pair>::iterator begin,
      std::vector<uint64_t_pair>::iterator end, std::uint8_t max_bits,
      T compare);  //  unary comparison function

  template <typename T>
  static std::vector<std::uint64_t> LongestSubsequence(  // only Match
      std::vector<uint64_t_pair>::const_iterator begin,
      std::vector<uint64_t_pair>::const_iterator end,
      T compare);  // binary comparison function

  std::uint32_t k_;
  std::uint32_t w_;
  std::uint32_t occurrence_;
  std::uint32_t m_;
  std::uint64_t g_;
  std::uint8_t n_;
  std::uint32_t best_n_;
  std::uint32_t reduce_win_sz_;
  bool robust_winnowing_;
  bool hpc_;
  std::vector<std::vector<minimizer_t>> minimizers_;
  std::vector<std::unordered_map<  // kmer -> (begin, count)
      uint128_t, std::pair<std::uint32_t, std::uint32_t>>>
      index_;
  std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
};

}  // namespace ram

#endif  // RAM_MINIMIZER_ENGINE_HPP_
