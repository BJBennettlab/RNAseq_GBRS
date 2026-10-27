#!/usr/bin/env python3
"""Download the GBRS reference files from their Zenodo record and verify them.

Files come from the record's own file list, each checked against the MD5 that
Zenodo publishes. Archives (bowtie index, transition probabilities) are unpacked
into the reference directory. The script writes reference_manifest.json with the
record, every downloaded file's MD5, and SHA-256 for every file now present.

An existing reference directory is reused only when every file still matches
its manifest entry; a mismatch stops the script instead of overwriting.

    python3 scripts/fetch_references.py REF_DIR [--record 15091310]
"""
import argparse
import hashlib
import json
import shutil
import tarfile
import urllib.request
from pathlib import Path

RECORD = "15091310"  # Reference Data for EMASE and GBRS, GRCm39 / Ensembl 105, 2025-03-26
NEEDED = [
    "bowtie_index.tar.gz",
    "transition_probabilities.tar.gz",
    "emase.fullTranscripts.info",
    "emase.pooled.fullTranscripts.info",
    "emase.gene2transcripts.tsv",
    "gbrs_emissions_all_tissues.avecs.npz",
    "ref.gene_pos.ordered_ensBuild_105.npz",
    "ref.genome_grid.GRCm39.tsv",
    "founder.hexcolor.info",
    "Mus_musculus.GRCm39.dna.primary_assembly.fa.fai",
]


def file_hash(path, name):
    h = hashlib.new(name)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def record_files(record):
    with urllib.request.urlopen(f"https://zenodo.org/api/records/{record}", timeout=120) as response:
        meta = json.load(response)
    files = {f["key"]: f for f in meta["files"]}
    missing = [n for n in NEEDED if n not in files]
    if missing:
        raise SystemExit(f"Record {record} lacks {missing}")
    return meta, files


def safe_extract(archive, dest):
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            target = (dest / member.name).resolve()
            if not str(target).startswith(str(dest.resolve()) + "/") or member.issym() or member.islnk():
                raise SystemExit(f"Unsafe archive member {member.name} in {archive.name}")
        tar.extractall(dest)


def link_fai(ref):
    """GBRS reconstruct reads chromosome lengths from GBRS_DATA/ref.fa.fai."""
    link = ref / "ref.fa.fai"
    if not link.exists():
        link.symlink_to("Mus_musculus.GRCm39.dna.primary_assembly.fa.fai")


def check_existing(ref):
    manifest = ref / "reference_manifest.json"
    data = json.loads(manifest.read_text())
    for rel, sha in data["sha256"].items():
        path = ref / rel
        if not path.is_file() or file_hash(path, "sha256") != sha:
            raise SystemExit(f"{path} differs from reference_manifest.json; use a new REF_DIR")
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ref_dir", type=Path)
    parser.add_argument("--record", default=RECORD)
    args = parser.parse_args()
    ref = args.ref_dir
    if (ref / "reference_manifest.json").exists():
        data = check_existing(ref)
        link_fai(ref)
        print(f"Reference directory verified: {len(data['sha256'])} files")
        return
    if ref.exists() and any(ref.iterdir()):
        raise SystemExit(f"{ref} is not empty and has no manifest; use a new directory")
    meta, files = record_files(args.record)
    ref.mkdir(parents=True, exist_ok=True)
    downloads = ref / ".downloads"
    downloads.mkdir(exist_ok=True)
    md5 = {}
    for name in NEEDED:
        entry = files[name]
        expected = entry["checksum"].removeprefix("md5:")
        part = downloads / (name + ".part")
        url = entry["links"]["self"]
        with urllib.request.urlopen(url, timeout=600) as response, part.open("wb") as out:
            shutil.copyfileobj(response, out, 1 << 22)
        got = file_hash(part, "md5")
        if got != expected:
            part.unlink()
            raise SystemExit(f"MD5 mismatch for {name}: {got} != {expected}")
        md5[name] = got
        if name.endswith(".tar.gz"):
            safe_extract(part, ref)
            part.rename(downloads / name)
        else:
            part.rename(ref / name)
        print(f"verified {name}")
    link_fai(ref)
    sha = {}
    for path in sorted(p for p in ref.rglob("*") if p.is_file() and not p.is_symlink() and ".downloads" not in p.parts):
        sha[str(path.relative_to(ref))] = file_hash(path, "sha256")
    manifest = {
        "record": args.record,
        "doi": meta.get("doi"),
        "title": meta["metadata"]["title"],
        "publication_date": meta["metadata"].get("publication_date"),
        "downloaded_md5": md5,
        "sha256": sha,
    }
    (ref / "reference_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    shutil.rmtree(downloads)
    print(f"Wrote reference_manifest.json: {len(sha)} files")


if __name__ == "__main__":
    main()
