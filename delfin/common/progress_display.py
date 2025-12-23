"""Thread-safe progress display with automatic TTY detection.

This module provides a simple, zero-dependency progress display that:
- Automatically detects TTY vs batch environments
- Uses carriage return for live updates in terminals
- Falls back to logging in batch/SLURM environments
- Thread-safe for use in background monitoring threads
"""

import sys
import threading
import logging
from typing import Optional

# Global registry to track active progress displays
_active_progress_displays = []
_registry_lock = threading.Lock()


class ProgressDisplay:
    """Thread-safe progress display with automatic TTY detection.

    Usage:
        progress = ProgressDisplay()
        progress.update("Processing: 50/100 (50%)")
        progress.finalize("Complete!")

    In TTY mode, updates appear on the same line.
    In batch mode, each update is on a new line.
    """

    def __init__(self, stream=None):
        """Initialize progress display.

        Args:
            stream: Output stream (default: sys.stderr)
        """
        self.stream = stream or sys.stderr
        self.is_tty = self.stream.isatty()
        self._lock = threading.Lock()
        self._last_line_length = 0
        self._active = True

        # Register this display in the global registry
        with _registry_lock:
            _active_progress_displays.append(self)

    def update(self, message: str, persistent: bool = False):
        """Update progress display.

        Args:
            message: Progress message to display
            persistent: If True, always write a newline (for important messages)
        """
        if not self._active:
            return

        with self._lock:
            try:
                if self.is_tty and not persistent:
                    # TTY mode: use carriage return for same-line update
                    # Clear previous line by padding with spaces
                    clear_spaces = ' ' * max(0, self._last_line_length - len(message))
                    self.stream.write(f'\r{message}{clear_spaces}')
                    self.stream.flush()
                    self._last_line_length = len(message)
                else:
                    # Batch mode or persistent message: write with newline
                    if self.is_tty and self._last_line_length > 0:
                        # Move to new line after progress display
                        self.stream.write('\n')
                    self.stream.write(f'{message}\n')
                    self.stream.flush()
                    self._last_line_length = 0
            except OSError:
                # Terminal went away (rare), disable further updates
                self._active = False

    def clear(self):
        """Clear the current progress line (TTY mode only)."""
        if not self._active:
            return

        with self._lock:
            if self.is_tty and self._last_line_length > 0:
                # Clear line by writing spaces and returning to start
                self.stream.write('\r' + ' ' * self._last_line_length + '\r')
                self.stream.flush()
                self._last_line_length = 0

    def finalize(self, final_message: Optional[str] = None, remove_from_registry: bool = True):
        """Finalize progress display (move to new line in TTY mode).

        Args:
            final_message: Optional final message to display
            remove_from_registry: If True, remove this display from the global registry
        """
        if not self._active:
            return

        with self._lock:
            if self.is_tty and self._last_line_length > 0:
                if final_message:
                    # Replace progress line with final message
                    clear_spaces = ' ' * max(0, self._last_line_length - len(final_message))
                    self.stream.write(f'\r{final_message}{clear_spaces}\n')
                else:
                    # Just move to next line
                    self.stream.write('\n')
                self.stream.flush()
                self._last_line_length = 0
            elif final_message:
                # Batch mode: write final message
                self.stream.write(f'{final_message}\n')
                self.stream.flush()

        # Remove from registry when fully finalized
        if remove_from_registry:
            with _registry_lock:
                try:
                    _active_progress_displays.remove(self)
                except ValueError:
                    pass  # Already removed


def clear_active_progress_lines():
    """Clear all active progress lines before logging output.

    This ensures that logger messages don't appear on the same line
    as progress displays in TTY mode.
    """
    with _registry_lock:
        for display in _active_progress_displays[:]:  # Copy list to avoid modification during iteration
            if display.is_tty and display._last_line_length > 0:
                # Finalize without removing from registry (so it can be resumed)
                display.finalize(remove_from_registry=False)


class ProgressAwareStreamHandler(logging.StreamHandler):
    """StreamHandler that clears progress displays before emitting log records.

    This prevents log messages from appearing on the same line as progress displays.
    """

    def emit(self, record):
        """Emit a log record, clearing active progress displays first."""
        try:
            # Clear any active progress lines before logging
            clear_active_progress_lines()
            # Emit the log record normally
            super().emit(record)
        except Exception:
            self.handleError(record)
