// Copyright (c) 2020 Robert Vaser

#include <getopt.h>

#include <bitset>
#include <cstdlib>
#include <iostream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/progress_bar.hpp"
#include "biosoup/timer.hpp"

#include "ram/minimizer_engine.hpp"

std::atomic<std::uint32_t> biosoup::Sequence::num_objects{0};

namespace {

const char* ram_version = RAM_VERSION;

static struct option options[] = {
    {"kmer-length", required_argument, nullptr, 'k'},
    {"window-length", required_argument, nullptr, 'w'},
    {"frequency-threshold", required_argument, nullptr, 'f'},
    {"Micromize", no_argument, nullptr, 'M'},
    {"m", required_argument, nullptr, 'm'},
    {"g", required_argument, nullptr, 'g'},
    {"n", required_argument, nullptr, 'n'},
    {"preset-options", required_argument, nullptr, 'x'},
    {"threads", required_argument, nullptr, 't'},
    {"version", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}};

std::unique_ptr<bioparser::Parser<biosoup::Sequence>> CreateParser(
    const std::string& path) {
  auto is_suffix = [](const std::string& s, const std::string& suff) {
    return s.size() < suff.size()
               ? false
               : s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
  };

  if (is_suffix(path, ".fasta") || is_suffix(path, ".fa") ||
      is_suffix(path, ".fasta.gz") || is_suffix(path, ".fa.gz")) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<
          bioparser::FastaParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }
  if (is_suffix(path, ".fastq") || is_suffix(path, ".fq") ||
      is_suffix(path, ".fastq.gz") || is_suffix(path, ".fq.gz")) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<
          bioparser::FastqParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  std::cerr << "[ram::CreateParser] error: file " << path
            << " has unsupported format extension (valid extensions: .fasta, "
            << ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)"
            << std::endl;
  return nullptr;
}

void Help() {
  // clang-format off
  std::cout
      << "usage: ram [options ...] <target> [<sequences>]\n"
         "\n"
         "  # default output is stdout\n"
         "  <target>/<sequences> \n"
         "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
         "\n"
         "  options will be applied sequentially as specified, example:\n"
         "  $ ram -w10 -k19 -w5 reads.fastq\n"
         "  will result in w = 5\n"
         "\n"
         "  options:\n"
         "    -k, --kmer-length <int>\n"
         "      default: 15\n"
         "      length of minimizers\n"
         "    -w, --window-length <int>\n"
         "      default: 5\n"
         "      length of sliding window from which minimizers are found\n"
         "    -f, --frequency-threshold <float>\n"
         "      default: 0.001\n"
         "      threshold for ignoring most frequent minimizers\n"
         "    -M, --Micromize\n"
         "      use only a portion of all minimizers\n"
         "    -m <int>\n"
         "      default: 100\n"
         "      discard chains with chaining score less than <int>\n"
         "    -g <int>\n"
         "      default: 10000\n"
         "      stop chain elongation if there are no minimizer withing <int>-BP\n"
         "    -n <int>\n"
         "      default: 4\n"
         "      discard chains consisting of less then <int> minimizers\n"
         "    -x, --preset-options ava|pb\n"
         "      default: none\n"
         "      preset options; applies multiple options at the same time;\n"
         "      this options will be overwritten if used with other options;\n"
         "      available preset options strings:\n"
         "          ava: all-vs-all alignment (-k19 -w5 -m100 -g10000 -n4)\n"
         "          map: read to reference mapping (-k19 -w10 -m40 -g5000 -n3)\n"
         "    -t, --threads <int>\n"
         "      default: 1\n"
         "      number of threads\n"
         "    --version\n"
         "      prints the version number\n"
         "    -h, --help\n"
         "      prints the usage\n";
  // clang-format on
}

}  // namespace

int main(int argc, char** argv) {
  std::uint32_t k = 15;
  std::uint32_t w = 5;
  double frequency = 0.001;
  bool micromize = false;
  std::uint32_t m = 100;
  std::uint64_t g = 10000;
  std::uint8_t n = 4;
  std::string preset = "";
  std::uint32_t num_threads = 1;

  std::vector<std::string> input_paths;

  const char* optstr = "k:w:f:Mm:g:n:x:t:h";
  char arg;
  // clang-format off
  while ((arg = getopt_long(argc, argv, optstr, options, nullptr)) != -1) {
    switch (arg) {
      case 'k': k = std::atoi(optarg); break;
      case 'w': w = std::atoi(optarg); break;
      case 'f': frequency = std::atof(optarg); break;
      case 'M': micromize = true; break;
      case 'm': m = std::atoi(optarg); break;
      case 'g': g = std::atoll(optarg); break;
      case 'n': n = std::atoi(optarg); break;
      case 'x':
        preset = optarg;
        if (preset == "ava") {
          k = 19, w = 5, m = 100, g = 10000, n = 4;
          break;
        } else if (preset == "map") {
          k = 19, w = 10, m = 40, g = 5000, n = 3;
          break;
        }
        Help();
        return 1;
      case 't': num_threads = std::atoi(optarg); break;
      case 'v': std::cout << ram_version << std::endl; return 0;
      case 'h': Help(); return 0;
      default: return 1;
    }
  }
  // clang-format on

  if (argc == 1) {
    Help();
    return 0;
  }

  std::cerr << "[ram::] using options: "
            << "k = " << k << ", w = " << w << ", f = " << frequency
            << ", M = " << micromize << ", m = " << m << ", g = " << g
            << ", n = " << (int)n << ", x = " << preset
            << ", t = " << num_threads << std::endl;

  for (auto i = optind; i < argc; ++i) {
    input_paths.emplace_back(argv[i]);
  }

  if (input_paths.empty()) {
    std::cerr << "[ram::] error: missing target file" << std::endl;
    return 1;
  }

  auto tparser = CreateParser(input_paths[0]);
  if (tparser == nullptr) {
    return 1;
  }

  bool is_ava = false;
  std::unique_ptr<bioparser::Parser<biosoup::Sequence>> sparser = nullptr;
  if (input_paths.size() > 1) {
    sparser = CreateParser(input_paths[1]);
    if (sparser == nullptr) {
      return 1;
    }
    is_ava = input_paths[0] == input_paths[1];
  } else {
    sparser = CreateParser(input_paths[0]);
    is_ava = true;
  }

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads);
  ram::MinimizerEngine minimizer_engine{k, w, m, g, n, thread_pool};

  biosoup::Timer timer{};

  while (true) {
    timer.Start();

    std::vector<std::unique_ptr<biosoup::Sequence>> targets;
    try {
      targets = tparser->Parse(1U << 30);
    } catch (std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return 1;
    }

    if (targets.empty()) {
      break;
    }

    std::cerr << "[ram::] parsed " << targets.size() << " targets "
              << std::fixed << timer.Stop() << "s" << std::endl;

    timer.Start();

    minimizer_engine.Minimize(targets.begin(), targets.end());
    minimizer_engine.Filter(frequency);

    std::cerr << "[ram::] minimized targets " << std::fixed << timer.Stop()
              << "s" << std::endl;

    std::uint64_t num_targets = biosoup::Sequence::num_objects;
    biosoup::Sequence::num_objects = 0;

    while (true) {
      timer.Start();

      std::vector<std::unique_ptr<biosoup::Sequence>> sequences;
      try {
        sequences = sparser->Parse(1U << 30);
      } catch (std::invalid_argument& exception) {
        std::cerr << exception.what() << std::endl;
        return 1;
      }

      if (sequences.empty()) {
        break;
      }

      std::vector<std::future<std::vector<biosoup::Overlap>>> futures;
      for (const auto& it : sequences) {
        futures.emplace_back(thread_pool->Submit(
            [&](const std::unique_ptr<biosoup::Sequence>& sequence)
                -> std::vector<biosoup::Overlap> {
              return minimizer_engine.Map(sequence, is_ava, is_ava, micromize);
            },
            std::ref(it)));
      }

      biosoup::ProgressBar bar{static_cast<std::uint32_t>(sequences.size()),
                               16};

      for (auto& it : futures) {
        std::uint64_t rhs_offset = targets.front()->id;
        std::uint64_t lhs_offset = sequences.front()->id;
        for (const auto& jt : it.get()) {
          // clang-format off
          std::cout << sequences[(jt.lhs_id >> 1) - lhs_offset]->name << "\t"
                    << sequences[(jt.lhs_id >> 1) - lhs_offset]->data.size() << "\t"  // NOLINT
                    << jt.lhs_begin << "\t"
                    << jt.lhs_end << "\t"
                    << (jt.lhs_id & 1 ? "+" : "-") << "\t"
                    << targets[jt.rhs_id - rhs_offset]->name << "\t"
                    << targets[jt.rhs_id - rhs_offset]->data.size() << "\t"
                    << jt.rhs_begin << "\t"
                    << jt.rhs_end << "\t"
                    << jt.score << "\t"
                    << std::max(
                          jt.lhs_end - jt.lhs_begin,
                          jt.rhs_end - jt.rhs_begin) << "\t"
                    << 255
                    << std::endl;
          // clang-format on
        }

        if (++bar) {
          std::cerr << "[ram::] mapped " << bar.event_counter() << " sequences "
                    << "[" << bar << "] " << std::fixed << timer.Lap() << "s"
                    << "\r";
        }
      }
      std::cerr << std::endl;
      timer.Stop();

      if (is_ava && biosoup::Sequence::num_objects == num_targets) {
        break;
      }
    }

    sparser->Reset();
    biosoup::Sequence::num_objects = num_targets;
  }

  std::cerr << "[ram::] " << timer.elapsed_time() << "s" << std::endl;

  return 0;
}
