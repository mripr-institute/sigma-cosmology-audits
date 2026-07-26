#!/usr/bin/env python3
"""Download and cryptographically verify this audit's public input catalogues."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parent
DEFAULT_MANIFEST = ROOT / "data_sources.json"
CHUNK = 8 * 1024 * 1024


def digest(path: Path, algorithm: str) -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as stream:
        while block := stream.read(CHUNK):
            value.update(block)
    return value.hexdigest()


def expected_hash(source: dict[str, object]) -> tuple[str, str]:
    for algorithm in ("sha256", "md5"):
        if source.get(algorithm):
            return algorithm, str(source[algorithm]).lower()
    raise ValueError(f"{source.get('id')}: no checksum is pinned")


def verify(path: Path, source: dict[str, object]) -> tuple[bool, str, str]:
    algorithm, expected = expected_hash(source)
    if not path.is_file():
        return False, algorithm, "missing"
    actual = digest(path, algorithm)
    return actual == expected, algorithm, actual


def download(source: dict[str, object], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        str(source["url"]), headers={"User-Agent": "MRIPR-data-custody/1.0"}
    )
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            prefix=f".{destination.name}.", suffix=".part",
            dir=destination.parent, delete=False,
        ) as output:
            temporary = Path(output.name)
            with urllib.request.urlopen(request, timeout=60) as response:
                shutil.copyfileobj(response, output, length=CHUNK)
            output.flush()
            os.fsync(output.fileno())
        valid, algorithm, actual = verify(temporary, source)
        if not valid:
            raise RuntimeError(
                f"checksum mismatch for {source['id']}: {algorithm}={actual}"
            )
        os.replace(temporary, destination)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def safe_extract(source: dict[str, object], archive: Path) -> None:
    extract_to = source.get("extract_to")
    if not extract_to:
        return
    target = ROOT / str(extract_to)
    suffix = str(source.get("member_suffix", ""))
    with zipfile.ZipFile(archive) as bundle:
        members = [item for item in bundle.infolist() if item.filename.endswith(suffix)]
        expected = int(source.get("expected_members", len(members)))
        if len(members) != expected:
            raise RuntimeError(
                f"{source['id']}: expected {expected} archive members, found {len(members)}"
            )
        target.mkdir(parents=True, exist_ok=True)
        for member in members:
            name = PurePosixPath(member.filename).name
            if not name or name in {".", ".."}:
                raise RuntimeError(f"unsafe archive member: {member.filename!r}")
            destination = target / name
            temporary = destination.with_name(f".{destination.name}.part")
            try:
                with bundle.open(member) as source_stream, temporary.open("wb") as output:
                    shutil.copyfileobj(source_stream, output, length=CHUNK)
                    output.flush()
                    os.fsync(output.fileno())
                if destination.exists():
                    if digest(destination, "sha256") != digest(temporary, "sha256"):
                        raise RuntimeError(
                            f"refused to overwrite mismatching extracted member: {destination}"
                        )
                    temporary.unlink()
                    continue
                os.replace(temporary, destination)
            finally:
                temporary.unlink(missing_ok=True)


def load_manifest(path: Path) -> list[dict[str, object]]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema_version") != 1 or not isinstance(document.get("sources"), list):
        raise ValueError("unsupported data-source manifest")
    return document["sources"]


def verify_extracted(source: dict[str, object]) -> tuple[bool, str]:
    expected = source.get("extracted_aggregate_sha256")
    extract_to = source.get("extract_to")
    if not expected or not extract_to:
        return False, "not configured"
    target = ROOT / str(extract_to)
    suffix = str(source.get("member_suffix", ""))
    paths = sorted(target.glob(f"*{suffix}"), key=lambda path: path.name)
    expected_members = int(source.get("expected_members", len(paths)))
    if len(paths) != expected_members:
        return False, f"expected {expected_members} members, found {len(paths)}"
    aggregate = hashlib.sha256()
    for path in paths:
        aggregate.update(path.name.encode("utf-8"))
        aggregate.update(b"\0")
        aggregate.update(digest(path, "sha256").encode("ascii"))
        aggregate.update(b"\0")
    actual = aggregate.hexdigest()
    return actual == str(expected).lower(), actual


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--source", action="append", help="operate only on this source id")
    args = parser.parse_args()
    selected = set(args.source or [])
    failures = 0
    for source in load_manifest(args.manifest):
        if selected and str(source["id"]) not in selected:
            continue
        destination = ROOT / str(source["path"])
        valid, algorithm, actual = verify(destination, source)
        extracted_valid, extracted_actual = verify_extracted(source)
        if valid:
            print(f"OK       {source['id']}  {algorithm}={actual}")
        elif args.verify_only and extracted_valid:
            print(
                f"OK       {source['id']}  extracted_aggregate_sha256={extracted_actual} "
                "(source archive not retained)"
            )
        elif args.verify_only:
            print(f"INVALID  {source['id']}  {algorithm}={actual}", file=sys.stderr)
            failures += 1
            continue
        else:
            if destination.exists():
                print(
                    f"REFUSED  {source['id']}: existing file has {algorithm}={actual}; "
                    "move it aside explicitly before downloading",
                    file=sys.stderr,
                )
                failures += 1
                continue
            print(f"DOWNLOAD {source['id']}  {source['url']}")
            try:
                download(source, destination)
            except Exception as error:
                print(f"FAILED   {source['id']}: {error}", file=sys.stderr)
                failures += 1
                continue
            print(f"OK       {source['id']}  verified and installed")
        if destination.is_file():
            try:
                safe_extract(source, destination)
            except Exception as error:
                print(f"FAILED   {source['id']} extraction: {error}", file=sys.stderr)
                failures += 1
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
