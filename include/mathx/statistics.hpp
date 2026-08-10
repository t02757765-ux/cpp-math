#ifndef MATHX_STATISTICS_HPP
#define MATHX_STATISTICS_HPP

#include "export.hpp"
#include <vector>

namespace mathx {

//! Arithmetic mean.
MATHX_API double mean(const std::vector<double>& x);
//! Sample variance (n-1 denominator).
MATHX_API double variance(const std::vector<double>& x);
//! Population variance (n denominator).
MATHX_API double variance_pop(const std::vector<double>& x);
//! Sample standard deviation.
MATHX_API double stddev(const std::vector<double>& x);
//! Population standard deviation.
MATHX_API double stddev_pop(const std::vector<double>& x);
//! Median.
MATHX_API double median(std::vector<double> x);
//! Minimum value.
MATHX_API double min(const std::vector<double>& x);
//! Maximum value.
MATHX_API double max(const std::vector<double>& x);
//! Sum.
MATHX_API double sum(const std::vector<double>& x);
//! Product.
MATHX_API double product(const std::vector<double>& x);
//! Geometric mean.
MATHX_API double geometric_mean(const std::vector<double>& x);
//! Harmonic mean.
MATHX_API double harmonic_mean(const std::vector<double>& x);
//! k-th percentile (linear interpolation between closest ranks).
MATHX_API double percentile(std::vector<double> x, double p);
//! Sample covariance of two arrays.
MATHX_API double covariance(const std::vector<double>& x,
                            const std::vector<double>& y);
//! Pearson correlation coefficient.
MATHX_API double correlation(const std::vector<double>& x,
                             const std::vector<double>& y);
//! Skewness of a sample.
MATHX_API double skewness(const std::vector<double>& x);
//! Excess kurtosis of a sample.
MATHX_API double kurtosis(const std::vector<double>& x);
//! Range = max - min.
MATHX_API double range(const std::vector<double>& x);

} // namespace mathx

#endif
