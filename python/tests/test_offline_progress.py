"""Tests for offline-driver progress reporting."""

from __future__ import annotations

from io import StringIO

from utahlsm_offline import _ProgressReporter


class _TTYBuffer(StringIO):
    """String buffer that presents itself as an interactive terminal."""

    def isatty(self) -> bool:
        """Reports interactive-terminal behavior."""
        return True


def test_interactive_progress_updates_one_terminal_line() -> None:
    """Interactive progress uses carriage returns and finishes with a newline."""
    stream = _TTYBuffer()
    timestamps = iter([0.0, 0.0, 5.0, 10.0])
    progress = _ProgressReporter(
        10,
        stream=stream,
        min_interval=0.0,
        clock=lambda: next(timestamps),
    )

    progress.update(0)
    progress.update(5)
    progress.update(10)

    rendered = stream.getvalue()
    assert rendered.count('\r') == 3
    assert '[**************--------------] ...Running:  50.0% (5/10)' in rendered
    assert 'ETA 00:05' in rendered
    assert rendered.endswith('\n')


def test_redirected_progress_emits_quarterly_milestones() -> None:
    """Redirected output remains readable and limited to coarse milestones."""
    stream = StringIO()
    current_time = 0.0

    def clock() -> float:
        return current_time

    progress = _ProgressReporter(20, stream=stream, clock=clock)
    for completed in range(21):
        current_time = float(completed)
        progress.update(completed)

    lines = stream.getvalue().splitlines()
    assert len(lines) == 5
    assert [line.split('%')[0].split()[-1] for line in lines] == [
        '0.0',
        '25.0',
        '50.0',
        '75.0',
        '100.0',
    ]


def test_progress_can_be_disabled() -> None:
    """Batch callers can suppress all progress output."""
    stream = _TTYBuffer()
    progress = _ProgressReporter(10, enabled=False, stream=stream)

    progress.update(10)

    assert stream.getvalue() == ''
