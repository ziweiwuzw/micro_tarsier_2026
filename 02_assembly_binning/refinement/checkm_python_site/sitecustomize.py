"""Raise Python's recursion limit for CheckM/DendroPy trees with many bins."""

from __future__ import absolute_import

import os
import sys


def _configured_limit():
    raw = os.environ.get("CHECKM_RECURSION_LIMIT", "100000")
    try:
        value = int(raw)
    except ValueError:
        return 100000
    return max(value, 1000)


sys.setrecursionlimit(_configured_limit())
