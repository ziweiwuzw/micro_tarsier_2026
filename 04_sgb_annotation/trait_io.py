#!/usr/bin/env python3
"""Small, dependency-free helpers for atomic trait-table exports."""

from __future__ import annotations

import csv
import gzip
import os
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def dict_writer(path: Path, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(temporary, "wt", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                    lineterminator="\n")
            writer.writeheader()
            yield writer
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def read_tsv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")
