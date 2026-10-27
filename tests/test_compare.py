import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("compare_gbrs_runs", ROOT / "scripts" / "compare_gbrs_runs.py")
compare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare)


def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    body = "locus\tA\tB\tC\tD\tE\tF\tG\tH\ttotal\tnotes\n"
    for g, a, b, note in rows:
        body += f"{g}\t{a}\t{b}\t0\t0\t0\t0\t0\t0\t{a + b}\t{note}\n"
    path.write_text(body)


class Compare(unittest.TestCase):
    def test_identical_and_changed_samples(self):
        with tempfile.TemporaryDirectory(prefix="gbrs_compare_codetestfile_", dir="/tmp") as d:
            root = Path(d)
            rows = [(f"g{i}", 10 * i + 5, 3 * i, "AB" if i % 2 else "AA") for i in range(1, 40)]
            for s in ("S1", "S2"):
                write(root / "new" / s / "either" / f"{s}.diploid.genes.expected_read_counts", rows)
            write(root / "old" / "S1.tsv", rows)
            write(root / "old" / "S2.tsv", [(g, b, a, "BB") for g, a, b, _ in rows])
            sheet = root / "samples.tsv"
            sheet.write_text("sample\nS1\nS2\nS3\n")
            out = root / "cmp" / "either_vs_old.tsv"
            compare.main(["--sheet", str(sheet), "--a", str(root / "new/{sample}/either/{sample}.diploid.genes.expected_read_counts"),
                          "--b", str(root / "old/{sample}.tsv"), "--out", str(out)])
            summary = json.loads(out.with_suffix(".json").read_text())
            self.assertEqual(summary["identical"], ["S1"])
            self.assertEqual(summary["missing"], {"a": ["S3"], "b": ["S3"]})
            lines = out.read_text().splitlines()
            self.assertTrue(lines[1].startswith("S1\t39\t39\t"))
            self.assertTrue(lines[2].endswith("\t0.0"))


if __name__ == "__main__":
    unittest.main()
