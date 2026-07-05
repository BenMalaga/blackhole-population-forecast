"""Stage-1 acquisition accounting: counts + bytes from manifests and the data tree.

Structure-only: prints file/row/byte counts, never sample values (firewall-safe).
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"


def tree_bytes(p: Path) -> int:
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file()) if p.exists() else 0


def main() -> None:
    log_p = DATA / "manifests" / "o4a_pe_acquisition_log.json"
    if log_p.exists():
        log = json.loads(log_p.read_text())
        evs = log["events"]
        n_pq = sum(len(e["extractions"]) for e in evs)
        n_samp = sum(x["n_samples"] for e in evs for x in e["extractions"])
        in_cut = sum(1 for e in evs if e.get("in_lvk_bbh_cut"))
        print(f"o4a-pe: {len(evs)}/86 events | {in_cut} in BBH cut | "
              f"{n_pq} parquet subsets | {n_samp} posterior samples (count only)")
    sel_p = DATA / "manifests" / "selection_extraction_manifest.json"
    if sel_p.exists():
        for m in json.loads(sel_p.read_text()):
            print(f"selection: {m['tag']}: {m['n_rows']} rows x {m['n_columns_kept']} cols "
                  f"({m['parquet_bytes']/1e6:.1f} MB parquet)")
    g_p = DATA / "manifests" / "gwtc3_tarball_stream_manifest.json"
    if g_p.exists():
        g = json.loads(g_p.read_text())
        print(f"gwtc3-plp: {len(g['extracted'])} members kept of "
              f"{g['members_seen_count']} seen; {g['compressed_bytes_read']/1e9:.2f} GB "
              f"compressed read; early_stop={g['early_stop']}")
        for e in g["extracted"]:
            print(f"   {e['member']} ({e['size']/1e6:.1f} MB) sha256={e['sha256'][:16]}…")
    for sub in sorted(p for p in DATA.iterdir() if p.is_dir()):
        print(f"data/{sub.name}: {tree_bytes(sub)/1e6:.1f} MB")
    print(f"TOTAL data/: {tree_bytes(DATA)/1e9:.3f} GB | .venv: "
          f"{tree_bytes(ROOT/'.venv')/1e9:.3f} GB")


if __name__ == "__main__":
    main()
