# Conditional locus statistics

These methods compare rates within one supplied DNA-copy model, tree, root
distribution, observation set, and event catalogue. They have not received
independent empirical calibration or a completed ablation study.

## Nested rate comparison

`locus-statistics --null-rate-group GROUP` fits a full model and a nested
model with the named free rate group fixed at zero. Other free groups are
re-estimated. When both fits resolve, the statistic is
`2 * (logLik_full - logLik_null)`, subject to the declared numerical
tolerance. Failed or unresolved fits yield no comparison statistic.
Multistart optimization and local curvature diagnostics describe numerical
behavior; they cannot establish a global maximum or scientific
identifiability. Profiles fix one rate at supplied values and refit nuisance
rates.

This comparison requires the same DNA units, tree, root, tip evidence,
observation coding, and opportunity catalogue. It cannot rank mapping
methods or state codings by raw likelihood when those choices alter observed
units or sampled loci.

## Fixed-catalogue bootstrap

With `--bootstrap-replicates B --seed N --sampling-design fixed_catalogue`,
the command simulates DNA observations under the fitted null and refits both
models. It keeps the supplied catalogue, tree, root, tip missingness, and
detection design fixed. Every requested replicate has a recorded status. A
failed replicate prevents a bootstrap P value from being reported. For
complete replicates, the one-sided Monte Carlo estimate uses
`(1 + number(T_sim >= T_observed)) / (B + 1)`.

This is plug-in calibration conditional on a fixed, prequalified catalogue.
It does not correct for selecting loci or event opportunities from the same
tip observations, and it does not establish Type-I error under an independent
biological truth. Tip observations in this interface are DNA calls only.

## What remains to be measured

An independent generator or externally specified truth is needed to estimate
Type-I error, power, ancestral-state accuracy, localization error, and
probability calibration. Discovery and missingness designs must be fixed
before simulation. The [ablation protocol](ablation_study.md) separates those
questions from internal numerical correctness checks. Selected empirical
examples cannot on their own calibrate error rates or P values.
