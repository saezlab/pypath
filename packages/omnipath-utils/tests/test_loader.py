"""Tests for the on-demand background loading module."""

import time

from omnipath_utils.db._loader import (
    _lock,
    _pending,
    is_pending,
    request_table,
)


class TestBackgroundLoader:
    def _clear_pending(self):
        with _lock:
            _pending.clear()

    def test_request_table_queues(self):
        """request_table returns True and marks the key as pending."""
        self._clear_pending()

        # Use a fake DB URL so the background thread will fail quickly
        result = request_table(
            'test_src',
            'test_tgt',
            9999,
            db_url='postgresql://fake/db',
        )
        assert result is True
        assert is_pending('test_src', 'test_tgt', 9999)

        # Let the thread fail and clean up
        time.sleep(1)
        self._clear_pending()

    def test_duplicate_request_is_noop(self):
        """Second request for same table returns False."""
        self._clear_pending()

        request_table(
            'dup_src',
            'dup_tgt',
            9999,
            db_url='postgresql://fake/db',
        )
        result2 = request_table(
            'dup_src',
            'dup_tgt',
            9999,
            db_url='postgresql://fake/db',
        )
        assert result2 is False

        time.sleep(1)
        self._clear_pending()

    def test_not_pending_initially(self):
        """Unknown table is not pending."""
        assert not is_pending('nonexistent', 'type', 9606)

    def test_pending_cleared_after_thread_finishes(self):
        """After the background thread completes (or fails), pending is cleared."""
        self._clear_pending()

        request_table(
            'clear_src',
            'clear_tgt',
            9999,
            db_url='postgresql://fake/db',
        )
        assert is_pending('clear_src', 'clear_tgt', 9999)

        # Wait for thread to finish (it will fail due to fake DB)
        time.sleep(3)

        assert not is_pending('clear_src', 'clear_tgt', 9999)
