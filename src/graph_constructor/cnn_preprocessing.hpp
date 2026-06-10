#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

namespace raven {

enum CnnChannel : std::size_t {
  kMultMin = 0,
  kMultMax = 1,
  kMultMean = 2,
  kMultMedian = 3,

  kAvgQualMin = 4,
  kAvgQualMax = 5,
  kAvgQualMean = 6,
  kAvgQualMedian = 7,

  kMinQualMin = 8,
  kMinQualMax = 9,
  kMinQualMean = 10,
  kMinQualMedian = 11,

  kNumChannels = 12
};

struct WindowStats {
  float min = 0.0f;
  float max = 0.0f;
  float mean = 0.0f;
  float median = 0.0f;
};

struct WindowedCnnInput {
  std::vector<float> X;   // channel-major: [12 * num_windows]
  std::int64_t num_channels = 12;
  std::int64_t num_windows = 0;
  bool defined = false;
};

template <typename T>
WindowStats ComputeWindowStats(
    const std::vector<T>& values,
    std::size_t begin,
    std::size_t end) {
  WindowStats stats{};

  if (begin >= end || end > values.size()) {
    return stats;
  }

  const std::size_t n = end - begin;
  T min_val = values[begin];
  T max_val = values[begin];
  double sum = 0.0;

  for (std::size_t i = begin; i < end; ++i) {
    const T v = values[i];
    if (v < min_val) min_val = v;
    if (v > max_val) max_val = v;
    sum += static_cast<double>(v);
  }

  stats.min = static_cast<float>(min_val);
  stats.max = static_cast<float>(max_val);
  stats.mean = static_cast<float>(sum / static_cast<double>(n));

  std::vector<T> tmp(values.begin() + begin, values.begin() + end);
  const std::size_t mid = n / 2;

  std::nth_element(tmp.begin(), tmp.begin() + mid, tmp.end());
  const T upper = tmp[mid];

  if ((n % 2) == 1) {
    stats.median = static_cast<float>(upper);
  } else {
    std::nth_element(tmp.begin(), tmp.begin() + (mid - 1), tmp.begin() + mid);
    const T lower = tmp[mid - 1];
    stats.median = static_cast<float>(
        (static_cast<double>(lower) + static_cast<double>(upper)) / 2.0);
  }

  return stats;
}

WindowedCnnInput ComputeWindowFeaturesInference(
    const std::vector<float>& mult,
    const std::vector<std::uint32_t>& avg_qual,
    const std::vector<std::uint32_t>& min_qual,
    std::size_t window_size = 10);

}  // namespace ram