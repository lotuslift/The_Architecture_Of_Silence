#!/usr/bin/env python3
"""Verify the portable repository checksum ledger."""

from __future__ import annotations

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "SHA256SUMS"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    failures: list[str] = []
    checked = 0
    for line_number, raw in enumerate(LEDGER.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        try:
            expected, relative = line.split(maxsplit=1)
        except ValueError:
            failures.append(f"line {line_number}: malformed checksum entry")
            continue
        relative = relative.lstrip("* ")
        path = ROOT / relative
        if not path.is_file():
            failures.append(f"missing: {relative}")
            continue
        actual = sha256(path)
        checked += 1
        if actual != expected:
            failures.append(f"mismatch: {relative}\n  expected {expected}\n  actual   {actual}")

    if failures:
        print("SHA-256 verification failed:")
        for failure in failures:
            print(f"- {failure}")
        return 1

    print(f"PASS: {checked} repository files match SHA256SUMS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
