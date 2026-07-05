"""Idempotent, disk-budgeted data acquisition for the out-of-sample BH population test.

All sources are open-access Zenodo records (CC-BY-4.0, anonymous download), verified
live 2026-06-10 via the Zenodo REST API. Every download is md5-verified against the
record's own file listing (snapshotted into ``data/zenodo_listings/`` and committed).

HARD DISK BUDGET (Stage-1 phase rule, enforced IN CODE)
    At most ``PHASE_BUDGET_BYTES`` (3.0 GB) of new files may exist at any moment
    (project ``data/`` tree + ``.venv``), and at least ``FREE_FLOOR_BYTES`` of disk
    must remain free after any download. Every download calls :func:`budget_check`
    first and the script refuses (exit, nothing fetched) if the action would exceed
    either limit. Discipline: download -> extract needed columns -> DELETE raw.
    Artifacts that cannot fit even transiently (single file > budget) are NOT
    downloaded; their staged plan is documented in data/README.md.

INTEGRITY FIREWALL (amended 2026-06-10, Stage-1 directive, documented, not silent)
    Original draft rule: no GWTC-4/5 per-event posterior *downloaded or opened*
    before PRE_REGISTRATION.md locks. Amended rule, per the Stage-1 acquisition
    directive (recorded here and in PRE_REGISTRATION_DRAFT.md / data/README.md):

    * Hyperposterior + event-posterior files are INPUT data and may be acquired
      pre-lock, BUT only through the BLIND mechanical extraction in this script:
      column subsets are written to parquet, raw HDF5 deleted; no statistic of the
      sample values is computed, printed, or plotted. The acquisition is logged
      (UTC timestamps, md5s, sha256s, columns) in data/manifests/.
    * NO out-of-sample predictive score, summary statistic, histogram, or plot of
      any GWTC-4/5 posterior may be produced before PRE_REGISTRATION.md is locked
      and committed. Scoring code must (and will) re-check that lock.
    * Tier ``gwtc5-pe`` additionally stays hard-guarded: its largest file (6.59 GB)
      exceeds the Stage-1 budget, and the lock should exist before O4b is touched.

Usage (run inside the project venv: ``.venv/bin/python -m src.fetch_data ...``):
    python -m src.fetch_data --plan               # show tiers + budget, fetch nothing
    python -m src.fetch_data --tier metadata      # ~25 MB: listings, event lists, docs
    python -m src.fetch_data --tier selection     # O4a injections + O3-era mixture ->
                                                  #   parquet, raw deleted (1.7 GB peak)
    python -m src.fetch_data --tier gwtc3-plp     # STREAM the 10.07 GB monolith, keep
                                                  #   only analyses/PowerLawPeak/* (~MBs)
    python -m src.fetch_data --tier o4a-pe        # BLIND per-event acquisition, one
                                                  #   file at a time (<0.6 GB transient)
    python -m src.fetch_data --tier gwtc4-pop     # REFUSED in Stage 1 (7.47 GB tar),
                                                  #   prints the staged plan
    python -m src.fetch_data --tier gwtc5         # REFUSED in Stage 1 (2.94 GB), plan
    python -m src.fetch_data --tier gwtc5-pe      # HARD-GUARDED (pre-reg lock + budget)
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import re
import shutil
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
VENV = ROOT / ".venv"
LISTINGS = DATA / "zenodo_listings"
MANIFESTS = DATA / "manifests"
SCRATCH = DATA / "pe_scratch"
PREREG = ROOT / "PRE_REGISTRATION.md"
PREREG_DRAFT = ROOT / "PRE_REGISTRATION_DRAFT.md"

ZENODO_FILE = "https://zenodo.org/records/{rec}/files/{name}?download=1"
ZENODO_API = "https://zenodo.org/api/records/{rec}"
GWOSC_API = "https://gwosc.org/eventapi/json/{catalog}/"

# ---- Stage-1 phase rules, enforced in code -------------------------------------------
PHASE_BUDGET_BYTES = 3_000_000_000   # hard cap on new files (data/ + .venv) at any moment
FREE_FLOOR_BYTES = 1_200_000_000     # never leave the machine with less free than this
CHUNK = 1 << 20

ALL_RECORDS = ["5655785", "16053484", "16740117", "16740128",
               "16911563", "19500064", "20292639"]

# Small pinned files (record, filename, md5, bytes, local subdir), md5s from the
# 2026-06-10 live verification. Anything not listed here is verified at download time
# against the snapshotted Zenodo listing JSON.
PINNED_METADATA = [
    # GWTC-4.0 population release: pins the exact FAR<1/yr BBH event-selection cut.
    ("16911563", "o4a_event_list.tar", "35ed4499b4abf8deccd53668a335bbf5", 10240, "gwtc4_pop"),
    ("16911563", "README.md", "2f991c0d9ab6e9fa3af155b7ce88d10b", 14274, "gwtc4_pop"),
    ("16911563", "figure_scripts.tar", None, 112640, "gwtc4_pop"),
    # GWTC-5.0 population release: the O4b-era selection cut.
    ("20292639", "Event_list.tar.gz", "a4d3a1054cccfca5d72b58f5673dcf0b", 1595, "gwtc5_pop"),
    ("20292639", "README.md", "2f3245f1c5cf48be80d16fc64bdb8114", 2499, "gwtc5_pop"),
    # Sensitivity-release notes: analysis time T_obs + injection conventions.
    ("16740117", "gwtc-4_o4a_sensitivity-estimates.md",
     "8ee18697e947c880764f78caa986e55d", 15148, "sensitivity"),
    ("19500064", "gwtc-5_o4ab_sensitivity-estimates.md",
     "457b02eca6088eb8cc23a97f37b01f8b", 16876, "sensitivity"),
    # GWTC-3 population release docs + small auxiliaries (the tutorial pins the
    # internal tarball paths of the frozen-fit hyperposterior).
    ("5655785", "README.md", "02489cc1cb8254dec5a1bb564b6ac059", 5140, "gwtc3_pop"),
    ("5655785", "py_requirements.txt", "09e6dd7ea79582bda70f93fcba971926", 115, "gwtc3_pop"),
    ("5655785", "table_data.tar.gz", "bd6f9984cc762cc44ff86f0520ce56d1", 4917, "gwtc3_pop"),
    ("5655785", "PowerLaw_Peak_DipBreak_Spline_Tutorial.html",
     "5c9f33a461948beba75f773dab59a601", 10631837, "gwtc3_pop"),
    # GWTC-4.0 PE release docs (allowed pre-lock per PRE_REGISTRATION_DRAFT §3:
    # tutorial/docs only, NOT PESummaryTable.hdf5, which contains per-event outcome
    # summaries and stays untouched until after lock).
    ("16053484", "GWTC4p0_PE_data_release.ipynb", None, 2880777, "gwtc4_pe"),
    ("16053484", "md5sums.txt", None, 52934, "gwtc4_pe"),
]

GWOSC_CATALOGS = ["GWTC-3-confident", "GWTC-4.0", "GWTC-5.0"]  # pin 4.0, NOT 4.1

# ---- GWTC-3 frozen fit: streaming extraction targets ----------------------------------
GWTC3_REC = "5655785"
GWTC3_TARBALL = "GWTC-3-population-data.tar.gz"
GWTC3_TARBALL_MD5 = "3c51561b5d8624210685b179c7d1f6ca"
GWTC3_TARBALL_BYTES = 10_070_439_649
# NOTE: the release README says `analysis/`, the tutorial says `analyses/`, accept
# both spellings (and any nesting) rather than gamble a 10 GB stream on one of them.
PLP_PATTERN = re.compile(r"(analys[ei]s|analysis)/PowerLawPeak/", re.IGNORECASE)
# Known member basenames (pinned from the release tutorial HTML, fetched in `metadata`).
PLP_KNOWN = {
    "o1o2o3_mass_c_iid_mag_iid_tilt_powerlaw_redshift_result.json",
    "o1o2o3_mass_c_iid_mag_iid_tilt_powerlaw_redshift_mass_data.h5",
    "o1o2o3_mass_c_iid_mag_iid_tilt_powerlaw_redshift_redshift_data.h5",
    "o1o2o3_mass_c_iid_mag_iid_tilt_powerlaw_redshift_magnitude_data.h5",
    "o1o2o3_mass_c_iid_mag_iid_tilt_powerlaw_redshift_orientation_data.h5",
}
PLP_MEMBER_CAP = 2_200_000_000  # refuse any single member bigger than this

# ---- Selection function: extract-and-delete ------------------------------------------
SELECTION_FILES = [
    # Smallest first: under disk pressure the mixture can succeed while the big O4a
    # set waits for the free-space floor.
    # GWTC-3-era cumulative mixture (reproduces the O3-era selection for the mock-
    # calibration gate). Cartesian-spin variant per PRE_REGISTRATION_DRAFT §12 TODO;
    # the polar twin is the same size if the lock decides otherwise.
    ("16740128", "mixture-semi_o1_o2-real_o3-cartesian_spins_20250503134659UTC.hdf",
     "e5d860f3f45863a657401b64ec31b5e8", 257339480),
    # O4a-only injection set (forward-models the detected O4a population).
    ("16740117", "samples-rpo4a_v2_20250503133839UTC-1366933504-23846400.hdf",
     "1cf34f975ebbb3c8cce461ee26467cd6", 1442289200),
]
# Injection table fields worth keeping. Schema verified 2026-06-10 from the files'
# own compound dtype (see data/manifests/selection_extraction_manifest.json): a single
# structured dataset `events` with ~119 (O4a) / ~35 (mixture) named fields.
INJ_COL_PATTERN = re.compile(
    r"(?i)^(mass[12]_source|q|z|redshift|luminosity_distance|"
    r"dluminosity_distance_dredshift|spin[12][xyz]|"
    r"spin[12]_(magnitude|polar_angle|azimuthal_angle)|chi_eff|chi_p|"
    r"lnpdraw_.*|weights?|semianalytic_.*|.*_far|.*_p_astro)$")
# Keep log-pdf draw components + weights at float64; physical params/FAR/p_astro at f32.
INJ_F64_PATTERN = re.compile(r"(?i)^(lnpdraw_|weights?$|semianalytic_weights)")

# ---- Outcome epoch O4a: blind per-event acquisition -----------------------------------
O4A_PE_REC = "16053484"
PE_FILE_SUFFIX = "-combined_PEDataRelease.hdf5"
WANTED_PE_COLS = [
    "mass_1_source", "mass_2_source", "chirp_mass_source", "total_mass_source",
    "mass_ratio", "redshift", "luminosity_distance", "comoving_distance",
    "a_1", "a_2", "tilt_1", "tilt_2", "cos_tilt_1", "cos_tilt_2", "chi_eff", "chi_p",
]
PE_META_KEYS = {"history", "version"}

GWTC5_PE_RECS = ["20348005", "20348006"]


# =========================================================================== utilities

def utcnow() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _tree_bytes(path: Path) -> int:
    if not path.exists():
        return 0
    total = 0
    for p in path.rglob("*"):
        try:
            if p.is_file() and not p.is_symlink():
                total += p.stat().st_size
        except OSError:
            pass
    return total


_VENV_BYTES_CACHE: int | None = None


def new_file_bytes() -> int:
    """Bytes of 'new files' this project currently holds (data tree + venv)."""
    global _VENV_BYTES_CACHE
    if _VENV_BYTES_CACHE is None:
        _VENV_BYTES_CACHE = _tree_bytes(VENV)
    return _tree_bytes(DATA) + _VENV_BYTES_CACHE


def budget_check(incoming_bytes: int, what: str, wait_s: int = 0) -> None:
    """Refuse any action that would break the Stage-1 disk rules. Exits loudly.

    A budget (held new files) violation is structural -> exit immediately. A free-disk
    floor violation can be transient (another workload owns the disk) -> optionally wait
    up to ``wait_s`` seconds for space to come back before giving up.
    """
    held = new_file_bytes()
    if held + incoming_bytes > PHASE_BUDGET_BYTES:
        sys.exit(
            f"[budget] REFUSED: {what} (+{incoming_bytes/1e9:.2f} GB) would put new-file "
            f"total at {(held + incoming_bytes)/1e9:.2f} GB > hard budget "
            f"{PHASE_BUDGET_BYTES/1e9:.1f} GB (currently held: {held/1e9:.2f} GB). "
            f"Delete extracted scratch or defer this artifact (see data/README.md)."
        )
    waited = 0
    while shutil.disk_usage(ROOT).free - incoming_bytes < FREE_FLOOR_BYTES:
        if waited >= wait_s:
            sys.exit(
                f"[budget] REFUSED: {what} (+{incoming_bytes/1e9:.2f} GB) would leave "
                f"{(shutil.disk_usage(ROOT).free - incoming_bytes)/1e9:.2f} GB free "
                f"(< floor {FREE_FLOOR_BYTES/1e9:.1f} GB)"
                + (f" after waiting {waited}s" if waited else "") + "."
            )
        if waited == 0:
            print(f"[wait] low disk for {what}: free "
                  f"{shutil.disk_usage(ROOT).free/1e9:.2f} GB; polling up to {wait_s}s",
                  flush=True)
        time.sleep(30)
        waited += 30


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=1, sort_keys=True))


def fetch_listing(rec: str, refresh: bool = False) -> dict:
    """Snapshot the Zenodo record listing (file ids, sizes, md5s) to a committed JSON."""
    LISTINGS.mkdir(parents=True, exist_ok=True)
    dest = LISTINGS / f"{rec}.json"
    if dest.exists() and not refresh:
        return json.loads(dest.read_text())
    with urllib.request.urlopen(ZENODO_API.format(rec=rec), timeout=60) as resp:
        record = json.load(resp)
    snap = {
        "record": rec,
        "fetched_utc": utcnow(),
        "title": record["metadata"]["title"],
        "publication_date": record["metadata"]["publication_date"],
        "doi": record.get("doi"),
        "files": [
            {"file_id": f["id"], "key": f["key"], "size": f["size"],
             "checksum": f["checksum"]}
            for f in sorted(record["files"], key=lambda x: x["key"])
        ],
    }
    write_json(dest, snap)
    print(f"[list] {rec}: {len(snap['files'])} files, "
          f"{sum(f['size'] for f in snap['files'])/1e9:.2f} GB -> {dest.relative_to(ROOT)}")
    return snap


def listing_entry(rec: str, name: str) -> dict:
    for f in fetch_listing(rec)["files"]:
        if f["key"] == name:
            return f
    sys.exit(f"[fail] {name} not found in Zenodo record {rec} listing")


def _download(rec: str, name: str, dest: Path, md5: str | None, size: int | None,
              download_wait_s: int = 0) -> None:
    """Download one Zenodo file -> dest, md5-verified, budget-checked, idempotent."""
    entry = listing_entry(rec, name)
    api_md5 = entry["checksum"].removeprefix("md5:")
    api_size = entry["size"]
    if md5 is not None and md5 != api_md5:
        sys.exit(f"[fail] {name}: pinned md5 {md5} != live listing md5 {api_md5}")
    if size is not None and size != api_size:
        sys.exit(f"[fail] {name}: pinned size {size} != live listing size {api_size}")
    md5, size = api_md5, api_size

    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        if dest.stat().st_size == size and _md5(dest) == md5:
            print(f"[skip] {dest.relative_to(ROOT)} (already verified)")
            return
        print(f"[redo] {dest.relative_to(ROOT)} size/md5 mismatch, refetching")
        dest.unlink()
    budget_check(size, f"download {name}", wait_s=download_wait_s)
    tmp = dest.with_suffix(dest.suffix + ".part")
    url = ZENODO_FILE.format(rec=rec, name=name)
    print(f"[get ] {name} ({size/1e6:.1f} MB) -> {dest.relative_to(ROOT)}")
    t0 = time.time()
    with urllib.request.urlopen(url, timeout=120) as resp, tmp.open("wb") as out:
        done = 0
        while True:
            chunk = resp.read(CHUNK)
            if not chunk:
                break
            out.write(chunk)
            done += len(chunk)
            if done % (1 << 26) < CHUNK:  # every ~64 MB: mid-download disk guard
                waited = 0
                while shutil.disk_usage(ROOT).free < 350_000_000 and waited < 300:
                    if waited == 0:
                        print(f"       [wait] disk critically low mid-download; pausing",
                              flush=True)
                    time.sleep(15)
                    waited += 15
                if shutil.disk_usage(ROOT).free < 350_000_000:
                    tmp.unlink()
                    sys.exit(f"[fail] {name}: aborted mid-download, disk critically low "
                             f"(partial file deleted; re-run to resume)")
            if size > 3e8 and done % (1 << 28) < CHUNK:
                print(f"       ... {done/1e9:.2f}/{size/1e9:.2f} GB "
                      f"({done/max(time.time()-t0, 1e-9)/1e6:.0f} MB/s)", flush=True)
    if tmp.stat().st_size != size:
        got = tmp.stat().st_size
        tmp.unlink()
        sys.exit(f"[fail] {name}: got {got} bytes, expected {size}")
    if _md5(tmp) != md5:
        tmp.unlink()
        sys.exit(f"[fail] {name}: md5 mismatch vs listing {md5}")
    tmp.rename(dest)


def _require_locked_prereg(tier: str) -> None:
    if not PREREG.exists():
        sys.exit(
            f"[firewall] tier '{tier}' is hard-guarded: {PREREG.name} does not exist.\n"
            f"Lock {PREREG_DRAFT.name} -> {PREREG.name}, git-commit it, then re-run."
        )
    text = PREREG.read_text(encoding="utf-8", errors="replace")
    if "NOT LOCKED" in text or "DRAFT" in text.splitlines()[0].upper():
        sys.exit(f"[firewall] {PREREG.name} is still marked DRAFT/NOT LOCKED, refusing '{tier}'.")


# ====================================================================== tier: metadata

def tier_metadata() -> None:
    for rec in ALL_RECORDS:
        fetch_listing(rec)
    for rec, name, md5, size, subdir in PINNED_METADATA:
        _download(rec, name, DATA / subdir / name, md5, size)
    out_dir = DATA / "gwosc"
    out_dir.mkdir(parents=True, exist_ok=True)
    for catalog in GWOSC_CATALOGS:
        dest = out_dir / f"{catalog}.json"
        if dest.exists():
            print(f"[skip] {dest.relative_to(ROOT)}")
            continue
        with urllib.request.urlopen(GWOSC_API.format(catalog=catalog), timeout=60) as resp:
            payload = json.load(resp)
        dest.write_text(json.dumps(payload, indent=1))
        print(f"[get ] GWOSC {catalog}: {len(payload['events'])} events "
              f"-> {dest.relative_to(ROOT)}")


# ============================================================ tier: gwtc3-plp (stream)

class _CountingReader:
    """Wrap an HTTP stream: count + md5 the compressed bytes, report progress."""

    def __init__(self, raw, total: int):
        self.raw, self.total = raw, total
        self.n = 0
        self.md5 = hashlib.md5()
        self.t0 = time.time()
        self._next_report = 1 << 30

    def read(self, size: int = -1) -> bytes:
        chunk = self.raw.read(size)
        if chunk:
            self.n += len(chunk)
            self.md5.update(chunk)
            if self.n >= self._next_report:
                self._next_report += 1 << 30
                rate = self.n / max(time.time() - self.t0, 1e-9) / 1e6
                print(f"       ... streamed {self.n/1e9:.2f}/{self.total/1e9:.2f} GB "
                      f"compressed ({rate:.0f} MB/s)", flush=True)
        return chunk


def tier_gwtc3_plp() -> None:
    """Stream the 10.07 GB GWTC-3 monolith; keep ONLY analyses/PowerLawPeak/* members.

    The raw tarball never touches disk (Stage-1 rule: artifact > budget => no full
    download). We decompress the HTTP stream on the fly, write whitelisted members,
    record a manifest of every member seen, and stop early once the PowerLawPeak
    directory has been fully passed. Early stop means the tarball's md5 cannot be
    fully verified; each extracted member is therefore pinned by sha256 here, and the
    extraction is reproducible by re-running this tier.
    """
    entry = listing_entry(GWTC3_REC, GWTC3_TARBALL)
    if entry["checksum"].removeprefix("md5:") != GWTC3_TARBALL_MD5:
        sys.exit("[fail] GWTC-3 tarball md5 changed upstream, re-verify before streaming")
    out_root = DATA / "gwtc3_pop"
    manifest_path = MANIFESTS / "gwtc3_tarball_stream_manifest.json"
    done_targets = {p.name for p in (out_root / "analyses" / "PowerLawPeak").glob("*")
                    if p.is_file()} if (out_root / "analyses" / "PowerLawPeak").exists() else set()
    if PLP_KNOWN <= done_targets and manifest_path.exists():
        print("[skip] PowerLawPeak members already extracted + manifest present")
        return

    url = ZENODO_FILE.format(rec=GWTC3_REC, name=GWTC3_TARBALL)
    print(f"[strm] {GWTC3_TARBALL} ({GWTC3_TARBALL_BYTES/1e9:.2f} GB compressed), "
          f"streaming, extracting analyses/PowerLawPeak/* only")
    members_seen: list[dict] = []
    extracted: list[dict] = []
    entered_plp = False
    early_stop = False
    with urllib.request.urlopen(url, timeout=120) as resp:
        counter = _CountingReader(resp, GWTC3_TARBALL_BYTES)
        tf = tarfile.open(fileobj=counter, mode="r|gz")
        for member in tf:
            members_seen.append({"name": member.name, "size": member.size})
            is_plp = bool(PLP_PATTERN.search(member.name)) and member.isfile()
            if is_plp:
                entered_plp = True
                if member.size > PLP_MEMBER_CAP:
                    print(f"[warn] {member.name} is {member.size/1e9:.2f} GB > cap, skipped")
                    continue
                budget_check(member.size, f"extract {member.name}", wait_s=240)
                rel = member.name[PLP_PATTERN.search(member.name).end():]
                dest = out_root / "analyses" / "PowerLawPeak" / rel
                dest.parent.mkdir(parents=True, exist_ok=True)
                src = tf.extractfile(member)
                h = hashlib.sha256()
                with dest.open("wb") as out:
                    for chunk in iter(lambda: src.read(CHUNK), b""):
                        h.update(chunk)
                        out.write(chunk)
                extracted.append({"member": member.name, "size": member.size,
                                  "sha256": h.hexdigest(),
                                  "path": str(dest.relative_to(ROOT))})
                print(f"[keep] {member.name} ({member.size/1e6:.1f} MB)")
            elif entered_plp:
                # tar members are grouped by directory; we've left PowerLawPeak.
                have = {Path(e["member"]).name for e in extracted}
                if PLP_KNOWN <= have:
                    early_stop = True
                    print(f"[strm] all PowerLawPeak targets captured after "
                          f"{counter.n/1e9:.2f} GB compressed, stopping early")
                    break
    write_json(manifest_path, {
        "tarball": GWTC3_TARBALL, "record": GWTC3_REC,
        "tarball_md5_pinned": GWTC3_TARBALL_MD5,
        "tarball_bytes": GWTC3_TARBALL_BYTES,
        "streamed_utc": utcnow(),
        "compressed_bytes_read": counter.n,
        "early_stop": early_stop,
        "compressed_md5_if_complete": counter.md5.hexdigest() if not early_stop else None,
        "members_seen_count": len(members_seen),
        "members_seen": members_seen,
        "extracted": extracted,
        "note": ("Raw tarball intentionally never stored (Stage-1 disk budget). "
                 "Extracted members pinned by sha256; manifest covers members seen "
                 "up to the early-stop point only."),
    })
    missing = PLP_KNOWN - {Path(e["member"]).name for e in extracted}
    if missing:
        sys.exit(f"[fail] streaming ended without capturing: {sorted(missing)}")
    print(f"[done] gwtc3-plp: {len(extracted)} members kept, manifest -> "
          f"{manifest_path.relative_to(ROOT)}")


# ================================================================ tier: selection (S)

def _find_injection_table(h5):
    """Find the injection table: the largest 1-D *structured* dataset in the file.

    Schema fact (verified 2026-06-10, recorded in the selection manifest): both the
    O4a set and the cumulative mixtures store one compound dataset named ``events``
    at the file root.
    """
    import h5py  # lazy: needs the project venv
    found = []
    def visit(name, obj):
        if (isinstance(obj, h5py.Dataset) and obj.ndim == 1
                and obj.dtype.names and obj.shape[0] > 1000):
            found.append((name, obj.shape[0], len(obj.dtype.names)))
    h5.visititems(visit)
    if not found:
        raise RuntimeError("no 1-D structured (compound-dtype) dataset found")
    return max(found, key=lambda t: t[1])[0]


def extract_injection_file(h5_path: Path, tag: str) -> dict:
    """Extract needed injection fields -> parquet (+attrs/structure JSON); delete raw.

    Raises (and the caller KEEPS the raw file) if nothing matches, an empty
    extraction must never silently destroy a 1.4 GB download.
    """
    import h5py
    import numpy as np
    import pandas as pd
    out_dir = DATA / "derived" / "selection"
    out_dir.mkdir(parents=True, exist_ok=True)
    info: dict = {"source_file": h5_path.name, "tag": tag, "extracted_utc": utcnow()}
    with h5py.File(h5_path, "r") as h5:
        table_name = _find_injection_table(h5)
        ds = h5[table_name]
        fields = list(ds.dtype.names)
        info["table_dataset"] = table_name
        info["n_fields_total"] = len(fields)
        info["all_fields"] = [{"name": f, "dtype": str(ds.dtype[f])} for f in fields]
        info["attrs_root"] = {k: _jsonable(v) for k, v in h5.attrs.items()}
        parent = table_name.rsplit("/", 1)[0]
        if parent:
            info["attrs_table_parent"] = {k: _jsonable(v)
                                          for k, v in h5[parent].attrs.items()}
        keep = [f for f in fields if INJ_COL_PATTERN.match(f)]
        if not keep:
            raise RuntimeError(f"no fields matched INJ_COL_PATTERN among {fields[:20]}...")
        n_rows = ds.shape[0]
        cols = {}
        for f in sorted(keep):
            arr = ds.fields(f)[...]  # partial I/O: reads this member only
            if arr.dtype.kind == "f" and not INJ_F64_PATTERN.match(f):
                arr = arr.astype(np.float32)
            if arr.dtype.kind in "fiub":
                cols[f] = arr
        if not cols:
            raise RuntimeError(f"matched fields {keep} but none were numeric")
        df = pd.DataFrame(cols)
    pq = out_dir / f"{tag}.parquet"
    df.to_parquet(pq, compression="zstd", index=False)
    info.update({"n_rows": int(n_rows), "columns_kept": sorted(cols),
                 "n_columns_kept": len(cols),
                 "parquet": str(pq.relative_to(ROOT)),
                 "parquet_bytes": pq.stat().st_size,
                 "parquet_sha256": _sha256(pq),
                 "source_md5": _md5(h5_path),
                 "source_bytes": h5_path.stat().st_size})
    h5_path.unlink()
    print(f"[xtr ] {h5_path.name}: {n_rows} rows x {len(cols)} cols -> "
          f"{pq.relative_to(ROOT)} ({pq.stat().st_size/1e6:.1f} MB); raw DELETED")
    return info


def _jsonable(v):
    import numpy as np
    if isinstance(v, (bytes, bytearray)):
        return v.decode(errors="replace")
    if isinstance(v, np.generic):
        return v.item()
    if isinstance(v, np.ndarray):
        return v.tolist() if v.size <= 64 else f"<array shape={v.shape} dtype={v.dtype}>"
    return v


def tier_selection() -> None:
    manifest_path = MANIFESTS / "selection_extraction_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else []
    done = {m["source_file"] for m in manifest}
    for rec, name, md5, size in SELECTION_FILES:
        if name in done:
            print(f"[skip] {name} (already extracted)")
            continue
        SCRATCH.mkdir(parents=True, exist_ok=True)
        raw = SCRATCH / name
        _download(rec, name, raw, md5, size, download_wait_s=3600)
        tag = ("o4a_injections" if "rpo4a" in name
               else "o1o2o3_mixture_cartesian" if "real_o3-cartesian" in name
               else Path(name).stem)
        info = extract_injection_file(raw, tag)
        info["zenodo_record"] = rec
        info["zenodo_file_id"] = listing_entry(rec, name)["file_id"]
        manifest.append(info)
        write_json(manifest_path, manifest)
    print(f"[done] selection: manifest -> {manifest_path.relative_to(ROOT)}")


# ============================================== tier: o4a-pe (BLIND acquisition, O4a)

def _o4a_pe_labels() -> dict[str, str]:
    """Per-event PE label used by the LVK GWTC-4.0 population analysis (BBH cut).

    Read from the pinned `o4a_event_list.tar` (selection metadata, allowed pre-lock):
    `events_list_bbh_only.txt` has lines `event,pe_label` (e.g. GW230605_065343,C00:Mixed).
    """
    src = DATA / "gwtc4_pop" / "o4a_event_list" / "events_list_bbh_only.txt"
    if not src.exists():
        sys.exit(f"[fail] {src.relative_to(ROOT)} missing, run --tier metadata, then "
                 f"untar o4a_event_list.tar (see data/README.md)")
    labels = {}
    for line in src.read_text().splitlines():
        if line.startswith("GW"):
            event, label = line.strip().split(",")
            labels[event] = label
    return labels


def extract_event_posterior(h5_path: Path, out_dir: Path, event: str,
                            lvk_label: str | None) -> dict:
    """BLIND mechanical extraction of one event's posterior column subsets.

    Extracts (a) the label the LVK GWTC-4.0 population analysis used for this event
    (from the pinned event list), and (b) the 'Mixed' analysis if different, so the
    pre-registration lock can choose either convention without re-downloading.
    Prints/records STRUCTURE only (group names, column names, row counts), never any
    statistic of the sample values. Caller deletes the raw HDF5 afterwards.
    """
    import h5py
    import numpy as np
    import pandas as pd
    with h5py.File(h5_path, "r") as h5:
        keys = list(h5.keys())
        cands = [k for k in keys if k.lower() not in PE_META_KEYS]
        wanted: list[str] = []
        if lvk_label is not None:
            exact = [k for k in cands if k == lvk_label]
            sub = [k for k in cands if lvk_label.split(":")[-1].lower() in k.lower()]
            if exact:
                wanted.append(exact[0])
            elif sub:
                wanted.append(sorted(sub)[0])
        mixed = sorted(k for k in cands if "mixed" in k.lower())
        if mixed and mixed[0] not in wanted:
            wanted.append(mixed[0])
        fallback_used = None
        if not wanted:
            # Events outside the BBH cut (e.g. NSBH) may carry no Mixed analysis at
            # all. They are pre-registered exclusions (PRE_REGISTRATION_DRAFT §7),
            # extract a best-effort single analysis for the exclusion table.
            pref = sorted(cands, key=lambda k: (
                0 if "xphm" in k.lower() else 1 if "seobnr" in k.lower() else 2, k))
            if pref:
                fallback_used = pref[0]
                wanted = [fallback_used]
        if not wanted:
            raise RuntimeError(f"no usable analysis group: lvk_label={lvk_label}, "
                               f"available={keys}")
        info = {
            "event": event,
            "source_file": h5_path.name,
            "analysis_groups_available": keys,
            "lvk_population_pe_label": lvk_label,
            "fallback_group_used": fallback_used,
            "extractions": [],
            "extracted_utc": utcnow(),
        }
        for group in wanted:
            g = h5[group]
            if "posterior_samples" not in g:
                raise RuntimeError(f"group '{group}' has no posterior_samples; "
                                   f"has {list(g.keys())}")
            table = g["posterior_samples"]
            if table.dtype.names:  # structured array layout
                names = list(table.dtype.names)
                cols_present = [c for c in WANTED_PE_COLS if c in names]
                if not cols_present:
                    raise RuntimeError(f"'{group}/posterior_samples' has none of the "
                                       f"wanted columns; fields: {names[:30]}")
                df = pd.DataFrame({c: table.fields(c)[...].astype(np.float32)
                                   for c in cols_present})
            else:  # 2-D matrix + parameter_names layout
                pn_key = next((k for k in g.keys() if "parameter_names" in k), None)
                if pn_key is None:
                    raise RuntimeError(f"'{group}/posterior_samples' is unstructured "
                                       f"and no parameter_names dataset exists")
                names = [n.decode() if isinstance(n, bytes) else str(n)
                         for n in g[pn_key][...]]
                cols_present = [c for c in WANTED_PE_COLS if c in names]
                if not cols_present:
                    raise RuntimeError(f"no wanted columns among {names[:30]}")
                mat = table[...]
                df = pd.DataFrame({c: np.asarray(mat[:, names.index(c)],
                                                 dtype=np.float32)
                                   for c in cols_present})
            safe_label = group.replace(":", "_").replace("/", "_")
            pq = out_dir / f"{event}__{safe_label}.parquet"
            pq.parent.mkdir(parents=True, exist_ok=True)
            df.to_parquet(pq, compression="zstd", index=False)
            info["extractions"].append({
                "group_used": group,
                "group_keys": list(g.keys()),
                "has_priors_group": "priors" in g,
                "n_samples": int(df.shape[0]),
                "columns_extracted": cols_present,
                "columns_wanted_missing": [c for c in WANTED_PE_COLS if c not in names],
                "parquet": str(pq.relative_to(ROOT)),
                "parquet_bytes": pq.stat().st_size,
                "parquet_sha256": _sha256(pq),
            })
    return info


def tier_o4a_pe(limit: int | None = None) -> None:
    """Acquire GWTC-4.0 per-event posterior SUBSETS, one file at a time, blind.

    Pre-lock acquisition is permitted under the amended Stage-1 firewall (see module
    docstring); each acquisition is logged. Files are processed in name order ==
    GPS-time order, so an interrupted/limited run is exactly the pre-registered
    'first N events by GPS time' subset.
    """
    locked = PREREG.exists()
    log_path = MANIFESTS / "o4a_pe_acquisition_log.json"
    log = json.loads(log_path.read_text()) if log_path.exists() else {
        "policy": ("BLIND PRE-LOCK ACQUISITION (amended firewall, 2026-06-10): column "
                   "subsets extracted mechanically; raw HDF5 deleted; no statistic of "
                   "any sample value computed, printed, or plotted; no predictive "
                   "score before PRE_REGISTRATION.md locks."),
        "record": O4A_PE_REC,
        "events": [],
    }
    if not locked:
        print("[firewall] PRE_REGISTRATION.md not locked: proceeding in BLIND "
              "acquisition mode (structure-only logging; scoring stays forbidden).")
    listing = fetch_listing(O4A_PE_REC)
    labels = _o4a_pe_labels()
    pe_files = [f for f in listing["files"] if f["key"].endswith(PE_FILE_SUFFIX)]
    pe_files.sort(key=lambda f: f["key"])  # name order == GPS-time order
    if limit is not None:
        pe_files = pe_files[:limit]
        print(f"[note] --limit {limit}: acquiring the first {limit} events by GPS time "
              f"(documented subset rule)")
    derived = DATA / "derived" / "o4a_pe"
    done = {e["source_file"] for e in log["events"]}
    SCRATCH.mkdir(parents=True, exist_ok=True)
    for i, f in enumerate(pe_files, 1):
        name = f["key"]
        if name in done:
            print(f"[skip] ({i}/{len(pe_files)}) {name}")
            continue
        # IGWN-GWTC4p0-<hash>_<n>-<EVENT>-combined_PEDataRelease.hdf5
        m = re.search(r"(GW\d{6}_\d{6})", name)
        if not m:
            sys.exit(f"[fail] cannot parse event name from '{name}'")
        event = m.group(1)
        raw = SCRATCH / name
        _download(O4A_PE_REC, name, raw, f["checksum"].removeprefix("md5:"), f["size"],
                  download_wait_s=1800)
        try:
            info = extract_event_posterior(raw, derived, event, labels.get(event))
        except Exception as exc:  # keep the failing file for mechanical inspection
            failed = SCRATCH / "FAILED" / name
            failed.parent.mkdir(parents=True, exist_ok=True)
            raw.rename(failed)
            sys.exit(f"[fail] extraction failed for {name}: {exc}\n"
                     f"       raw kept at {failed.relative_to(ROOT)} for structural "
                     f"inspection (budget: inspect then delete).")
        raw.unlink()
        info.update({"zenodo_file_id": f["file_id"],
                     "source_md5_verified": f["checksum"].removeprefix("md5:"),
                     "source_bytes": f["size"],
                     "in_lvk_bbh_cut": event in labels,
                     "prereg_locked_at_acquisition": locked})
        log["events"] = [e for e in log["events"] if e["source_file"] != name] + [info]
        write_json(log_path, log)
        groups = [x["group_used"] for x in info["extractions"]]
        ns = [x["n_samples"] for x in info["extractions"]]
        print(f"[blnd] ({i}/{len(pe_files)}) {event}: groups={groups}, "
              f"n_samples={ns}; raw DELETED")
    n_ev = len(log["events"])
    n_pq = sum(len(e["extractions"]) for e in log["events"])
    print(f"[done] o4a-pe: {n_ev} events, {n_pq} parquet subsets; log -> "
          f"{log_path.relative_to(ROOT)}")


# ================================================================ deferred tiers

def tier_gwtc4_pop() -> None:
    sys.exit(
        "[budget] REFUSED in Stage 1: analyses_BBH.tar is a 7.47 GB monolithic tar "
        "(Zenodo 16911563, md5 ba554251bacbda3979206ba673c644f5) > 3.0 GB phase budget.\n"
        "Staged plan (documented in data/README.md): once disk frees up, "
        "download -> untar only the PowerLawPeak/Default popsummary HDF5 -> delete tar; "
        "or stream-untar (plain tar, range-requestable) extracting selected members.\n"
        "Already fetched now (metadata tier): o4a_event_list.tar (the selection-cut "
        "pin), README.md, figure_scripts.tar, and the full file listing with md5s "
        "(data/zenodo_listings/16911563.json)."
    )


def tier_gwtc5() -> None:
    sys.exit(
        "[budget] DEFERRED in Stage 1: samples-rpo4ab-...-clipped.hdf is 2.94 GB "
        "(Zenodo 19500064, md5 a1106c27ec6cfd906231613523f7b174); it fits the budget "
        "alone but not alongside the Stage-1 working set. Staged plan: fetch with the "
        "same extract-columns-then-delete discipline (expected parquet ~0.6 GB) right "
        "before the O4b epoch work starts. The release notes "
        "(gwtc-5_o4ab_sensitivity-estimates.md) are already fetched by --tier metadata."
    )


def tier_gwtc5_pe() -> None:
    _require_locked_prereg("gwtc5-pe")
    sys.exit(
        "[budget] REFUSED in Stage 1 even though the pre-registration is locked: the "
        "largest GWTC-5.0 PE file (GW240615_113620, 6.59 GB; records 20348005/20348006, "
        "55.41 GB total) exceeds the 3.0 GB phase budget even one-at-a-time. Staged "
        "plan: when the disk reservation lifts, stream per-event with ~7 GB "
        "transient scratch using the same blind extract-and-delete loop as o4a-pe."
    )


# ========================================================================== plan/main

def print_plan() -> None:
    print(__doc__)
    held = new_file_bytes()
    free = shutil.disk_usage(ROOT).free
    print(f"Current new-file total (data/ + .venv): {held/1e9:.2f} GB of "
          f"{PHASE_BUDGET_BYTES/1e9:.1f} GB budget; disk free: {free/1e9:.2f} GB "
          f"(floor {FREE_FLOOR_BYTES/1e9:.1f} GB)")


def main() -> None:
    try:
        os.nice(19)
    except OSError:
        pass
    tiers = {
        "metadata": tier_metadata,
        "selection": tier_selection,
        "gwtc3-plp": tier_gwtc3_plp,
        "o4a-pe": None,  # handled below (takes --limit)
        "gwtc4-pop": tier_gwtc4_pop,
        "gwtc5": tier_gwtc5,
        "gwtc5-pe": tier_gwtc5_pe,
    }
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tier", choices=sorted(tiers))
    ap.add_argument("--plan", action="store_true", help="print the plan and exit")
    ap.add_argument("--limit", type=int, default=None,
                    help="o4a-pe: acquire only the first N events by GPS time")
    args = ap.parse_args()
    if args.plan or not args.tier:
        print_plan()
        return
    DATA.mkdir(exist_ok=True)
    if args.tier == "o4a-pe":
        tier_o4a_pe(limit=args.limit)
    else:
        tiers[args.tier]()
    held = new_file_bytes()
    print(f"[disk] new-file total now {held/1e9:.2f} GB / {PHASE_BUDGET_BYTES/1e9:.1f} GB "
          f"budget; free {shutil.disk_usage(ROOT).free/1e9:.2f} GB")
    print("[done]")


if __name__ == "__main__":
    main()
