# rquest 1.1.1

* Changed the default sample quantile type from 8 to 7 to match R's `quantile()`.

* When `dist` is omitted, `q.test()` and `qineq()` now use a lognormal working
  distribution for QOR estimation in samples with fewer than 50 non-missing
  observations. The rule applies separately to each sample. The lognormal
  option requires strictly positive observations.

* Clarified the documentation for arguments passed to `qcov()` and renamed
  objects in the `q.test()` examples for clarity.

# rquest 1.1.0

* Added `qden()` for quantile density estimation using QOR-based and density-based methods.

* Added `qor()` for computing quantile optimality ratios using the generalized lambda distribution, lognormal, normal, exponential, and user-supplied quantile functions.

* Added `qrcov()` for covariance matrix estimation for ratios of linear combinations of quantile estimators.

* Added `rcv.test()` for one- and two-sample inference for robust coefficients of variation based on the MAD or IQR.

* Extended `qcov()` to use the new quantile density estimation framework, including flexible distributional choices for QOR-based estimation.

* Extended `qineq()` to support additional quantile-based inequality measures and user-specified measures.

* Updated `q.test()` to improve quantile-based inference and argument handling.

* Improved input validation, documentation, and error and warning messages throughout the package.

* Expanded the automated test suite to cover the new and updated functions.

# rquest 1.0.5

* Changed the help files for u (qcov and q.test) and u2 (q.test)

# rquest 1.0.4

# rquest 1.0.3

# rquest 1.0.2

* Title and description changed according to CRAN requirements

# rquest 1.0.1

* Added a entries to `NEWS.md` file to track changes to the package.

# rquest 1.0.0

* Initial CRAN submission.
