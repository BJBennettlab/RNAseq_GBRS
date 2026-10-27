import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


sheet = load("sample_sheet")


class SampleSheet(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="gbrs_sheet_codetestfile_", dir="/tmp")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.reads = self.root / "reads"
        for s in ["S_1", "S_10", "S_100"]:
            d = self.reads / s
            d.mkdir(parents=True)
            for mate in ("R1", "R2"):
                (d / f"{s}_{mate}.fastq.gz").write_bytes(s.encode())

    def write_list(self, rows):
        path = self.root / "list.tsv"
        path.write_text("sample\tsex\tgeneration\n" + "".join("\t".join(r) + "\n" for r in rows))
        return path

    def test_prefix_names_keep_their_own_reads(self):
        rows = sheet.build(self.write_list([("S_1", "M", "43"), ("S_10", "F", "43"), ("S_100", "M", "42")]), self.reads)
        for sample, r1, r2, *_ in rows:
            self.assertEqual(Path(r1).name, f"{sample}_R1.fastq.gz")
            self.assertEqual(Path(r1).read_bytes(), sample.encode())

    def test_missing_reads_are_refused_not_matched_by_prefix(self):
        for p in (self.reads / "S_1").iterdir():
            p.unlink()
        with self.assertRaisesRegex(ValueError, "missing or empty"):
            sheet.build(self.write_list([("S_1", "M", "43"), ("S_10", "F", "43"), ("S_100", "M", "43")]), self.reads)

    def test_unlisted_and_duplicate_samples_are_refused(self):
        with self.assertRaisesRegex(ValueError, "not in the sample list"):
            sheet.build(self.write_list([("S_1", "M", "43")]), self.reads)
        self.assertEqual(len(sheet.build(self.write_list([("S_1", "M", "43")]), self.reads, allow_unlisted=True)), 1)
        with self.assertRaisesRegex(ValueError, "twice"):
            sheet.build(self.write_list([("S_1", "M", "43"), ("S_1", "M", "43")]), self.reads, allow_unlisted=True)

    def test_sex_and_generation_are_checked(self):
        with self.assertRaisesRegex(ValueError, "sex"):
            sheet.build(self.write_list([("S_1", "unknown", "43")]), self.reads, allow_unlisted=True)
        with self.assertRaisesRegex(ValueError, "generation"):
            sheet.build(self.write_list([("S_1", "F", "")]), self.reads, allow_unlisted=True)


class Combine(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="gbrs_combine_codetestfile_", dir="/tmp")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.samples = self.root / "samples"
        lines = ["sample\tr1\tr2\tr1_bytes\tr2_bytes\tsex\tgeneration"]
        for i, s in enumerate(["S_1", "S_10"]):
            d = self.samples / s / "either"
            d.mkdir(parents=True)
            body = "locus\tA\tB\tC\tD\tE\tF\tG\tH\ttotal\tnotes\n"
            body += f"g1\t{i}\t1\t0\t0\t0\t0\t0\t0\t{i + 1}\tAB\ng2\t0\t0\t2\t0\t0\t0\t0\t0\t2\tCC\n"
            out = d / f"{s}.diploid.genes.expected_read_counts"
            out.write_text(body)
            done = {"outputs": {f"either/{out.name}": hashlib.sha256(out.read_bytes()).hexdigest()}, "mate_rules": ["either"],
                    "bowtie": {"R1": {"reads": 10, "aligned": 9}, "R2": {"reads": 10, "aligned": 8}}, "versions": "test"}
            (d.parent / "done.json").write_text(json.dumps(done))
            lines.append(f"{s}\tr1\tr2\t1\t1\tM\t43")
        self.sheet = self.root / "samples.tsv"
        self.sheet.write_text("\n".join(lines) + "\n")

    def run_combine(self):
        return subprocess.run([sys.executable, str(ROOT / "scripts/combine_gbrs.py"), str(self.sheet),
                               str(self.samples), str(self.root / "combined"), "either"], capture_output=True, text=True)

    def test_tables_follow_the_sheet(self):
        result = self.run_combine()
        self.assertEqual(result.returncode, 0, result.stderr)
        counts = (self.root / "combined/gbrs_counts.tsv").read_text().splitlines()
        self.assertEqual(counts, ["gene_id\tS_1\tS_10", "g1\t1\t2", "g2\t2\t2"])
        self.assertIn("g1\tAB\tAB", (self.root / "combined/gbrs_diplotypes.tsv").read_text())

    def test_missing_rule_is_refused(self):
        result = subprocess.run([sys.executable, str(ROOT / "scripts/combine_gbrs.py"), str(self.sheet),
                                 str(self.samples), str(self.root / "combined_both"), "both"], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no results for the mate rule both", result.stderr)

    def test_unfinished_or_changed_sample_is_refused(self):
        (self.samples / "S_10" / "either" / "S_10.diploid.genes.expected_read_counts").write_text("changed")
        self.assertNotEqual(self.run_combine().returncode, 0)
        (self.samples / "S_10" / "done.json").unlink()
        result = self.run_combine()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not finished", result.stderr)


if __name__ == "__main__":
    unittest.main()
