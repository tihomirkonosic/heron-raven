#include "cnn_preprocessing.hpp"

#include <stdexcept>

namespace raven {

WindowedCnnInput ComputeWindowFeaturesInference(
    const std::vector<float>& mult,
    const std::vector<std::uint32_t>& avg_qual,
    const std::vector<std::uint32_t>& min_qual,
    std::size_t window_size) {
  if (window_size < 1) {
    throw std::invalid_argument("window_size must be >= 1");
  }

  const std::size_t n = mult.size();
  if (avg_qual.size() != n || min_qual.size() != n) {
    throw std::invalid_argument(
        "mult, avg_qual, and min_qual must have the same length");
  }

  WindowedCnnInput out{};

  if (n == 0) {
    out.defined = true;
    out.num_windows = 0;
    return out;
  }

  // Offset is fixed to 0 for inference.
  const std::size_t n_full = n / window_size;
  const std::size_t rem = n % window_size;
  const std::size_t n_windows = n_full + (rem ? 1 : 0);

  out.X.assign(kNumChannels * n_windows, 0.0f);
  out.num_windows = static_cast<std::int64_t>(n_windows);
  out.defined = true;

  auto write_x = [&](std::size_t channel, std::size_t w, float value) {
    out.X[channel * n_windows + w] = value;
  };

  for (std::size_t w = 0; w < n_full; ++w) {
    const std::size_t begin = w * window_size;
    const std::size_t end = begin + window_size;

    const WindowStats m = ComputeWindowStats(mult, begin, end);
    const WindowStats aq = ComputeWindowStats(avg_qual, begin, end);
    const WindowStats mq = ComputeWindowStats(min_qual, begin, end);

    write_x(kMultMin, w, m.min);
    write_x(kMultMax, w, m.max);
    write_x(kMultMean, w, m.mean);
    write_x(kMultMedian, w, m.median);

    write_x(kAvgQualMin, w, aq.min);
    write_x(kAvgQualMax, w, aq.max);
    write_x(kAvgQualMean, w, aq.mean);
    write_x(kAvgQualMedian, w, aq.median);

    write_x(kMinQualMin, w, mq.min);
    write_x(kMinQualMax, w, mq.max);
    write_x(kMinQualMean, w, mq.mean);
    write_x(kMinQualMedian, w, mq.median);
  }

  if (rem) {
    const std::size_t w = n_full;
    const std::size_t begin = n_full * window_size;
    const std::size_t end = n;

    const WindowStats m = ComputeWindowStats(mult, begin, end);
    const WindowStats aq = ComputeWindowStats(avg_qual, begin, end);
    const WindowStats mq = ComputeWindowStats(min_qual, begin, end);

    write_x(kMultMin, w, m.min);
    write_x(kMultMax, w, m.max);
    write_x(kMultMean, w, m.mean);
    write_x(kMultMedian, w, m.median);

    write_x(kAvgQualMin, w, aq.min);
    write_x(kAvgQualMax, w, aq.max);
    write_x(kAvgQualMean, w, aq.mean);
    write_x(kAvgQualMedian, w, aq.median);

    write_x(kMinQualMin, w, mq.min);
    write_x(kMinQualMax, w, mq.max);
    write_x(kMinQualMean, w, mq.mean);
    write_x(kMinQualMedian, w, mq.median);
  }

  return out;
}

}  // namespace ram