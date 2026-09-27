"""Throttled phase-local progress events for copy-pair staging."""
from __future__ import annotations

import time


_REPORT_EVERY = 128
_REPORT_SECONDS = 30.0


def print_progress_event(event):
    fields = [f"derive_progress_{event['event']}", f"phase={event['phase']}"]
    for name in (
        "processed_occurrence_pairs", "total_occurrence_pairs",
        "processed_copy_groups", "total_copy_groups",
        "processed_phase_items", "total_phase_items",
        "elapsed_seconds", "phase_eta_seconds",
        "phase_elapsed_seconds",
    ):
        if name in event:
            fields.append(f"{name}={event[name]}")
    print(" ".join(fields), flush=True)


class ProgressReporter:
    def __init__(self, callback, total_occurrence_pairs):
        self.callback = callback
        self.total_occurrence_pairs = int(total_occurrence_pairs)
        self.started = time.monotonic()
        self.phase = None
        self.phase_started = None
        self._last_report = 0
        self._last_report_time = None
        self._total_phase_items = None
        self._total_copy_groups = None

    def _event(self, event, **values):
        if self.callback is None:
            return
        payload = {
            "event": event,
            "phase": self.phase,
            "processed_occurrence_pairs": 0,
            "total_occurrence_pairs": self.total_occurrence_pairs,
            "elapsed_seconds": round(time.monotonic() - self.started, 3),
            **values,
        }
        self.callback(payload)

    def start(self, phase, *, total_phase_items=None, processed_occurrence_pairs=0,
              total_copy_groups=None):
        self.phase = str(phase)
        self.phase_started = time.monotonic()
        self._last_report = 0
        self._last_report_time = self.phase_started
        self._total_phase_items = total_phase_items
        self._total_copy_groups = total_copy_groups
        self._event(
            "phase_start",
            total_phase_items=total_phase_items,
            processed_occurrence_pairs=int(processed_occurrence_pairs),
            total_copy_groups=total_copy_groups,
            processed_copy_groups=0 if total_copy_groups is not None else None,
        )

    def update(self, processed, total, *, processed_occurrence_pairs=None,
               processed_copy_groups=None, force=False):
        processed, total = int(processed), int(total)
        if self.callback is None:
            return
        now = time.monotonic()
        if (
            not force and processed != total
            and processed - self._last_report < _REPORT_EVERY
            and self._last_report_time is not None
            and now - self._last_report_time < _REPORT_SECONDS
        ):
            return
        self._last_report = processed
        self._last_report_time = now
        phase_elapsed = now - self.phase_started
        values = {
            "processed_phase_items": processed,
            "total_phase_items": total,
            "phase_elapsed_seconds": round(phase_elapsed, 3),
        }
        if processed_occurrence_pairs is not None:
            values["processed_occurrence_pairs"] = int(processed_occurrence_pairs)
        if processed_copy_groups is not None:
            values["processed_copy_groups"] = int(processed_copy_groups)
        if processed > 0 and phase_elapsed > 0:
            values["phase_eta_seconds"] = round(max(0.0, phase_elapsed * (total - processed) / processed), 3)
            values["eta_scope"] = "current_phase"
        self._event("progress", **values)

    def complete(self, *, processed_occurrence_pairs=None, processed_copy_groups=None):
        values = {}
        if self._total_phase_items is not None:
            values["processed_phase_items"] = int(self._total_phase_items)
            values["total_phase_items"] = int(self._total_phase_items)
        if processed_occurrence_pairs is not None:
            values["processed_occurrence_pairs"] = int(processed_occurrence_pairs)
        if processed_copy_groups is not None:
            values["processed_copy_groups"] = int(processed_copy_groups)
        if self._total_copy_groups is not None:
            values["total_copy_groups"] = int(self._total_copy_groups)
            if processed_copy_groups is None:
                values["processed_copy_groups"] = int(self._total_copy_groups)
        self._event("phase_complete", **values)
