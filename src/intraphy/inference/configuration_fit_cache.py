"""Fit-local reuse of unit-rate sparse configuration generators."""
from __future__ import annotations

import numpy as np

from ..structure.edits import EDIT_KINDS
from .configuration_model import RateModel
from .configuration_sparse import sparse_generator
from .kernel_cache import KernelCache


class CommonRateKernels:
    """Lazily scale unit-rate templates for one exact state space."""

    def __init__(self, space, cache: KernelCache, cache_lock=None):
        self._space = space
        self._cache = cache
        self._cache_lock = cache_lock
        self._namespace = object()
        self._unit_model = RateModel({kind: 1. for kind in EDIT_KINDS})

    def generator(self, space, model, origins, child):
        if space is not self._space:
            raise ValueError("Sparse template cache is bound to a different state space")
        if not space.complete:
            raise ValueError("state_space_incomplete: refusing a renormalized truncated generator")
        mu = next(iter(model.rates.values()))
        if any(rate != mu for rate in model.rates.values()):
            raise ValueError("Sparse templates require equal edit rates")
        key = (self._namespace, tuple(sorted(
            material for material, origin in origins.items() if origin == child)))
        if self._cache_lock is None:
            template = self._cache.get_or_compute(
                key, lambda: sparse_generator(space, self._unit_model, origins, child))
        else:
            with self._cache_lock:
                template = self._cache.get_or_compute(
                    key, lambda: sparse_generator(space, self._unit_model, origins, child))
        coefficient = mu * model.scale * (
            model.foreground_multiplier if child in model.foreground else 1.)
        if not np.isfinite(coefficient) or coefficient < 0:
            raise ArithmeticError("Sparse generator contains an invalid transition rate")
        q = template.copy()
        q.data *= coefficient
        if not np.isfinite(q.data).all():
            raise ArithmeticError("Sparse generator contains nonfinite values")
        q.eliminate_zeros()
        return q
