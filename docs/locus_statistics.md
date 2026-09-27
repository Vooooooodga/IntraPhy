# Locus comparison statistics and calibration

These procedures define conditional maximum-likelihood comparisons for a fixed
locus catalogue and a parametric-bootstrap calibration interface. Type-I error,
power, and estimation accuracy have not been measured under independent
generative simulations.

## Nested comparison

`compare_locus_rates` compares a full rate model with a nested model in which
the named free rate groups are fixed to zero. Both fits use the same supplied
loci, tree, root prior, observed tip calls, state space, and event catalogue.
The likelihood-ratio statistic is `2 * (logLik_full - logLik_null)` when both
optimizations are usable. The comparison routine retains the raw likelihood
difference and its declared numerical tolerance; differences within that
tolerance are represented by statistic zero in either direction. Failed or
unresolved fits do not yield a statistic.

The optimizer uses multistart fits. The null fit fixes the tested rate groups
at zero and re-estimates the remaining free nuisance rates; the full fit also
uses the fitted null nuisance values as an additional starting point. The
caller-supplied start-scale grid defines initial-value multipliers only; it is
neither a parameter interval nor an upper bound. The zero-rate boundary is
retained in the null likelihood. An unsuccessful optimizer, or an unresolved
higher-likelihood solution, prevents production of a comparison statistic.

Profile calculations likewise fix a target rate at each requested value and
re-fit the remaining nuisance rates. The supplied profile grid is a set of
evaluation points, not an inferred interval or parameter bound. Boundary
likelihoods are retained; a profile point with an unsuccessful fit or an
unresolved higher-likelihood solution is reported without a likelihood-ratio
statistic.

## Bootstrap contract

`bootstrap_locus_comparison(bundle, null_rate_groups, replicates, seed, *,
sampling_design, ...)` simulates each replicate under the
fitted null rates, then refits null and full models. The caller must explicitly
declare `sampling_design="fixed_catalogue"`. This operation conditions on the
catalogue, tree, root prior, and observed material-missingness mask; it does
not resample loci, correct ascertainment, or support a catalogue selected by
the same observed results. DNA calls are simulated from the latent material
state using the supplied sensitivity/specificity values. Tip-level feature or
path observations are outside this DNA-only interface. Latent catalogue
features may remain in the state space, provided no feature/path observation
data are supplied.

Each replicate receives an independent NumPy
`SeedSequence([seed, replicate_index])`; this makes generated replicates
independent of worker scheduling. Replicate computations are process-parallel, while each
fit itself uses one worker. The statistic comparison already applies its
declared likelihood tolerance; bootstrap tail ties are counted with exact
`T_sim >= T_observed`, without adding a second ad hoc tolerance.

For `B` requested replicates, the reported one-sided Monte Carlo estimate is
`(1 + number(T_sim >= T_observed)) / (B + 1)`. Its smallest attainable value
and nominal resolution are `1/(B+1)`. `completed_replicates` means successful
replicates; the `replicates` field retains one record for every attempted
replicate, including failures. Any failed simulated dataset or fit makes
`p_value` unavailable (`None`), while preserving diagnostics. If the observed
comparison itself fails, no simulations are attempted and the result retains
the failed comparison object and requested replicate count.

This is a plug-in parametric-bootstrap calibration at estimated null rates,
conditional on the declared fixed design. It is not an exact finite-sample test
of a composite null, and no asymptotic chi-square P value is substituted.

## Calibration and validation still required

The implementation has not been assessed with an independent generative
truth. A future bounded calibration should include:

1. Null scenarios spanning boundary rates, interior nuisance-rate values,
   short and long trees, and sparse versus dense opportunities; compare the
   empirical rejection frequency with its Monte Carlo uncertainty.
2. Alternative scenarios with prespecified rate contrasts and event
   opportunities; quantify power and estimation behavior separately from
   type-I calibration.
3. Missingness scenarios with fixed masks of differing coverage and, in a
   separate design, masks generated independently of the latent states.
4. Detection-error scenarios spanning declared sensitivity/specificity
   values, including misspecified detection assumptions.
5. Model-misspecification scenarios (for example, rate heterogeneity or
   transition processes not represented by the fitted model).

The truth generator for these assessments must be independently implemented
or otherwise independently specified. Reusing this module's sampler as both
the sole truth generator and the method under evaluation would test internal
consistency, not calibration under independent generative truth.

## Method references

- NumPy documentation, “Parallel random number generation” (SeedSequence and
  independent parallel streams):
  <https://numpy.org/doc/stable/reference/random/parallel.html>
- SciPy documentation, `scipy.optimize.minimize` L-BFGS-B options:
  <https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html>
- Cohen, O., Ashkenazy, H., Belinky, F., Huchon, D. and Pupko, T. (2010),
  “GLOOME: gain loss mapping engine,” *Bioinformatics* 26(22): 2914–2915.
  <https://doi.org/10.1093/bioinformatics/btq549>
