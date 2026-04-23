from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable


def format_elapsed(seconds: float) -> str:
    total_seconds = max(0, int(round(seconds)))
    hours, remainder = divmod(total_seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours > 0:
        return f"{hours:02d}:{minutes:02d}:{secs:02d}"
    return f"{minutes:02d}:{secs:02d}"


@dataclass
class ProgressTracker:
    label: str
    total: int
    progress_callback: Callable[[str], None] | None
    start_time: float = field(default_factory=time.perf_counter)

    def emit(self, current: int, detail: str | None = None) -> None:
        if self.progress_callback is None or self.total <= 0:
            return

        current = max(0, min(int(current), int(self.total)))
        elapsed = time.perf_counter() - self.start_time
        fraction = current / float(self.total)
        pct = 100.0 * fraction
        eta = (elapsed / current) * (self.total - current) if current > 0 else float("nan")

        parts = [
            f"{self.label}: {current}/{self.total}",
            f"({pct:5.1f}%)",
            f"elapsed={format_elapsed(elapsed)}",
        ]
        if current > 0 and current < self.total:
            parts.append(f"eta={format_elapsed(eta)}")
        if detail:
            parts.append(detail)
        self.progress_callback(" | ".join(parts))
