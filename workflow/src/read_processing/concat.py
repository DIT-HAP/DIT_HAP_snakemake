"""Concatenation of per-sample insertion counts and annotations.

Holds the module documentation for ``concat_counts_and_annotations.py``. That
script's only function, ``concatenate``, reads every one of its inputs
(``counts_files``, ``annotation_files``, ``output_counts``,
``output_annotations``) straight from its ``Config`` dataclass and is pure
orchestration plus file I/O with no separable domain logic to extract (per
Phase 4 convention 2): loading each count/annotation file, concatenating them,
and writing the two combined tables. It dissolves entirely into that script's
``main()`` and contributes nothing here.

Input
-----
- Per-sample count and annotation DataFrames, as prepared by the invoking
  script.

Output
------
- None. This module intentionally holds no functions or dataclasses.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""
