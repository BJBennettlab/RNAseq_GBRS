import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("record_sample", ROOT / "scripts" / "record_sample.py")
record = importlib.util.module_from_spec(spec)
spec.loader.exec_module(record)

# Written by bowtie 1.3.1 in the GBRS 1.0.1 container (October 2026).
LOG = """# reads processed: 15717680
# reads with at least one alignment: 12161917 (77.38%)
# reads that failed to align: 3555763 (22.62%)
Reported 201696652 alignments
"""


class Record(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="gbrs_record_codetestfile_", dir="/tmp")
        self.addCleanup(self.temp.cleanup)
        self.d = Path(self.temp.name) / "S_1.partial"
        for rule in ("either", "both"):
            (self.d / rule).mkdir(parents=True)
            for name in record.NEEDED:
                (self.d / rule / f"S_1.{name}").write_text(f"{rule} {name}")
        for mate in ("R1", "R2"):
            (self.d / f"S_1.bowtie_{mate}.log").write_text(LOG)

    def run_record(self):
        return record.record(self.d, "S_1", "/r/S_1_R1.fastq.gz", "/r/S_1_R2.fastq.gz", "M", "43", "t.npz",
                             "image@sha256:x", "gbrs 1.0.0", "either both", "abc123")

    def test_bowtie_1_3_log(self):
        data = self.run_record()
        self.assertEqual(data["bowtie"]["R1"], {"reads": 15717680, "aligned": 12161917})
        saved = json.loads((self.d / "done.json").read_text())
        self.assertEqual(saved["pipeline_revision"], "abc123")
        self.assertIn("either/S_1.genotypes.tsv", saved["outputs"])
        self.assertEqual(saved["mate_rules"], ["either", "both"])

    def test_older_log_wording(self):
        text = LOG.replace("at least one alignment", "at least one reported alignment")
        self.assertEqual(record.bowtie_counts(text), {"reads": 15717680, "aligned": 12161917})

    def test_missing_output_or_unequal_mates_stop(self):
        (self.d / "both" / "S_1.genotypes.tsv").unlink()
        with self.assertRaisesRegex(SystemExit, "missing output both/S_1.genotypes.tsv"):
            self.run_record()
        self.assertFalse((self.d / "done.json").exists())
        (self.d / "both" / "S_1.genotypes.tsv").write_text("x")
        (self.d / "S_1.bowtie_R2.log").write_text(LOG.replace("15717680", "15717681"))
        with self.assertRaisesRegex(SystemExit, "read counts differ"):
            self.run_record()


if __name__ == "__main__":
    unittest.main()
