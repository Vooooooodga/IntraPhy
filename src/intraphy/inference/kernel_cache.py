"""Run-local bounded caches; eviction changes work, never the probability model."""
from __future__ import annotations
from collections import OrderedDict
import numpy as np


class KernelCache:
    def __init__(self, maximum_bytes=128*1024*1024):
        self.maximum_bytes = maximum_bytes
        self.entries = OrderedDict()
        self.bytes = 0
        self.hits = self.misses = 0

    @staticmethod
    def _size(value):
        seen = set()
        def visit(item):
            if isinstance(item, np.ndarray):
                if id(item) in seen:
                    return 0
                seen.add(id(item))
                return item.nbytes
            if all(hasattr(item, name) for name in ("data", "indices", "indptr")):
                if id(item) in seen:
                    return 0
                seen.add(id(item))
                return sum(visit(getattr(item, name))
                           for name in ("data", "indices", "indptr"))
            if isinstance(item, dict):
                return sum(visit(v) for v in item.values())
            if isinstance(item, (tuple, list)):
                return sum(visit(v) for v in item)
            return 0
        return visit(value)

    def get_or_compute(self, key, compute):
        if key in self.entries:
            self.hits += 1
            value, size = self.entries.pop(key)
            self.entries[key] = value, size
            return value
        self.misses += 1
        value = compute()
        size = self._size(value)
        if size <= self.maximum_bytes:
            while self.entries and self.bytes+size > self.maximum_bytes:
                _, (_, previous) = self.entries.popitem(last=False)
                self.bytes -= previous
            self.entries[key] = value, size
            self.bytes += size
        return value
