import csv
import gzip
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("validate_gbrs", ROOT / "scripts" / "validate_gbrs.py")
validate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validate)

GENES = [validate.XIST, *validate.Y_GENES] + [f"ENSMUSG{n:011d}" for n in range(1, 2201)]
CHROMS = ["1", "2", "X"]
POINTS = 40  # per chromosome


def mosaic(seed):
    """A random founder mosaic: one founder pair per block of eight grid points."""
    rng = np.random.default_rng(seed)
    out = np.zeros((POINTS * len(CHROMS), 8))
    for start in range(0, len(out), 8):
        a, b = rng.integers(0, 8, size=2)
        out[start:start + 8, a] += 0.5
        out[start:start + 8, b] += 0.5
    return out


def expression(seed, sex):
    rng = np.random.default_rng(1000 + seed)
    base = np.random.default_rng(7).gamma(2.0, 200.0, size=len(GENES))
    values = base * rng.lognormal(0, 0.6, size=len(GENES))
    values[0] = 3000.0 if sex == "F" else 5.0
    values[1:5] = 0.0 if sex == "F" else 120.0
    return values


class Validate(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="gbrs_validate_codetestfile_", dir="/tmp")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        grid = self.root / "grid.tsv"
        with grid.open("w") as f:
            f.write("marker\tchr\tpos\tcM\tbp\n")
            for c in CHROMS:
                for k in range(POINTS):
                    bp = 3_000_000 + k * 1_000_000
                    f.write(f"{c}_{bp}\t{c}\t{bp}\t0\t{bp}\n")
        self.grid = grid
        # Mice m1..m5. Library S3 holds mouse 1's reads (the 2025 failure); S5 is labeled female but is male.
        self.libraries = {"S1": (1, "M"), "S2": (2, "F"), "S3": (1, "M"), "S4": (4, "F"), "S5": (5, "M")}
        sheet_sex = {"S1": "M", "S2": "F", "S3": "M", "S4": "F", "S5": "F"}
        self.samples = self.root / "samples"
        lines = ["sample\tr1\tr2\tr1_bytes\tr2_bytes\tsex\tgeneration"]
        for s, (mouse, sex) in self.libraries.items():
            d = self.samples / s / "either"
            d.mkdir(parents=True)
            np.savetxt(d / f"{s}.gbrs.interpolated.genoprobs.tsv", mosaic(mouse), fmt="%.6f",
                       header="\t".join("ABCDEFGH"), delimiter="\t")
            with (d / f"{s}.diploid.genes.expected_read_counts").open("w") as f:
                f.write("locus\tA\tB\tC\tD\tE\tF\tG\tH\ttotal\tnotes\n")
                for g, v in zip(GENES, expression(mouse, sex)):
                    f.write(f"{g}\t{v}\t0\t0\t0\t0\t0\t0\t0\t{v}\tAA\n")
            (d.parent / "done.json").write_text(json.dumps({"mate_rules": ["either"], "pipeline_revision": "test"}))
            lines.append(f"{s}\tr1\tr2\t1\t1\t{sheet_sex[s]}\t43")
        self.sheet = self.root / "samples.tsv"
        self.sheet.write_text("\n".join(lines) + "\n")
        self.list = self.root / "list.tsv"
        self.list.write_text("sample\tsex\tgeneration\tmouse_id\n" + "".join(
            f"{s}\t{sheet_sex[s]}\t43\tm{s[1]}\n" for s in self.libraries))
        self.dna = self.root / "dna.tsv.gz"
        with gzip.open(self.dna, "wt") as f:
            f.write("mouse_id\tchr\tpos_mb\tA\tB\tC\tD\tE\tF\tG\tH\n")
            for m in range(1, 6):
                probs = mosaic(m)
                for i, c in enumerate(CHROMS):
                    for k in range(0, POINTS, 2):
                        row = probs[i * POINTS + k]
                        f.write(f"m{m}\t{c}\t{3 + k}\t" + "\t".join(f"{x:.3f}" for x in row) + "\n")
        self.star = self.root / "star.tsv"
        with self.star.open("w") as f:
            f.write("gene_id\t" + "\t".join(self.libraries) + "\n")
            values = {s: expression(m, sex) for s, (m, sex) in self.libraries.items()}
            for i, g in enumerate(GENES):
                f.write(f"{g}.5\t" + "\t".join(str(round(values[s][i] * 1.3)) for s in self.libraries) + "\n")

    def run_checks(self, *extra):
        out = self.root / "out"
        code = validate.main(["--sheet", str(self.sheet), "--samples", str(self.samples), "--rule", "either",
                              "--grid", str(self.grid), "--out", str(out), "--every", "1", *extra])
        return code, json.loads((out / "summary.json").read_text()), out

    def read(self, path):
        with open(path) as f:
            return {r["sample"]: r for r in csv.DictReader(f, delimiter="\t")}

    def test_two_labels_with_one_mouse_are_flagged(self):
        code, summary, _ = self.run_checks()
        self.assertEqual(code, 1)
        pairs = summary["checks"]["duplicates"]["flagged_pairs"]
        self.assertEqual([(p["a"], p["b"]) for p in pairs], [("S1", "S3")])
        self.assertEqual(summary["checks"]["dna"]["status"], "skipped")

    def test_genotype_group_marks_expected_matches(self):
        self.list.write_text("sample\tsex\tgeneration\tgenotype_group\n"
                             "S1\tM\t43\tstrainA\nS2\tF\t43\t\nS3\tM\t43\tstrainA\nS4\tF\t43\tstrainA\nS5\tF\t43\t\n")
        _, summary, _ = self.run_checks("--sample-list", str(self.list))
        pairs = {(p["a"], p["b"]) for p in summary["checks"]["duplicates"]["flagged_pairs"]}
        self.assertEqual(pairs, {("S1", "S4"), ("S3", "S4")})

    def test_sex_mismatch(self):
        _, summary, out = self.run_checks()
        self.assertEqual(summary["checks"]["sex"]["failed"],
                         [{"sample": "S5", "sheet_sex": "F", "expression_sex": "M"}])
        self.assertEqual(self.read(out / "sex.tsv")["S2"]["expression_sex"], "F")

    def test_dna_identity(self):
        _, summary, out = self.run_checks("--sample-list", str(self.list), "--dna", str(self.dna))
        rows = self.read(out / "dna_identity.tsv")
        self.assertEqual(rows["S1"]["best_mouse"], "m1")
        self.assertEqual(rows["S4"]["status"], "ok")
        self.assertEqual([f["sample"] for f in summary["checks"]["dna"]["failed"]], ["S3"])
        self.assertEqual(summary["checks"]["dna"]["failed"][0]["best_mouse"], "m1")

    def test_star_identity(self):
        _, summary, out = self.run_checks("--star", str(self.star))
        rows = self.read(out / "star_identity.tsv")
        self.assertEqual(rows["S2"]["status"], "ok")
        self.assertEqual(rows["S4"]["best_library"], "S4")
        self.assertEqual(rows["S3"]["best_library"] in ("S1", "S3"), True)

    def test_unfinished_sample_fails_the_run(self):
        (self.samples / "S4" / "done.json").unlink()
        code, summary, _ = self.run_checks()
        self.assertEqual(code, 1)
        self.assertEqual(summary["checks"]["finished"]["unfinished"], {"S4": "no done.json"})


if __name__ == "__main__":
    unittest.main()
