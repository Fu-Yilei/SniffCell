import stat
import tempfile
import unittest
from pathlib import Path

from sniffcell.discover.sv_post_processing import (
    _filter_sample_specific_by_ad,
    _resolve_args,
    sv_post_processing_main,
)


def _make_exec(path: Path, body: str) -> None:
    path.write_text(body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def _build_two_group_env(root: Path) -> tuple[Path, Path, dict[str, str]]:
    sample_dir = root / "samplex"
    split_dir = sample_dir / "deconv" / "deconv_requested_group_splits"
    split_dir.mkdir(parents=True)
    reference = root / "ref.fa"
    reference.write_text(">chr1\nA\n", encoding="utf-8")
    for group in ("A", "B"):
        bam = split_dir / f"{group}.bam"
        bai = split_dir / f"{group}.bam.bai"
        bam.write_text("", encoding="utf-8")
        bai.write_text("", encoding="utf-8")
        (split_dir / f"{group}.sniffles.vcf.gz").write_text("", encoding="utf-8")
        (split_dir / f"{group}.sniffles.vcf.gz.tbi").write_text("", encoding="utf-8")
    (split_dir / "requested_group_splits.tsv").write_text(
        "requested_group\tbam_path\tread_summary_path\n"
        f"A\t{split_dir / 'A.bam'}\t{split_dir / 'A.tsv'}\n"
        f"B\t{split_dir / 'B.bam'}\t{split_dir / 'B.tsv'}\n",
        encoding="utf-8",
    )

    bin_dir = root / "bin"
    bin_dir.mkdir()
    tools = {}
    for tool in ("bcftools", "truvari", "kanpig"):
        path = bin_dir / tool
        tools[tool] = str(path)
    _make_exec(
        Path(tools["bcftools"]),
        """#!/usr/bin/env python3
import sys
from pathlib import Path
args = sys.argv[1:]
sub = args[0]
if sub == "view":
    out = Path(args[args.index("-o") + 1]); out.parent.mkdir(parents=True, exist_ok=True); out.write_text("vcf\\n", encoding="utf-8")
elif sub == "index":
    target = Path(args[-1]); Path(str(target) + ".tbi").write_text("idx\\n", encoding="utf-8")
elif sub == "concat":
    out = Path(args[args.index("-o") + 1]); out.parent.mkdir(parents=True, exist_ok=True); out.write_text("vcf\\n", encoding="utf-8")
elif sub == "sort":
    out = Path(args[args.index("-o") + 1]); out.parent.mkdir(parents=True, exist_ok=True)
    src = Path(args[-1]); out.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
else:
    raise SystemExit(sub)
""",
    )
    _make_exec(
        Path(tools["truvari"]),
        """#!/usr/bin/env python3
import sys
from pathlib import Path
args = sys.argv[1:]
raw = Path(args[args.index("-o") + 1])
removed = Path(args[args.index("-c") + 1])
raw.parent.mkdir(parents=True, exist_ok=True)
raw.write_text("##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\nchr1\\t10\\tsv1\\tA\\t<DEL>\\t60\\tPASS\\tNumCollapsed=1\\n", encoding="utf-8")
removed.write_text("##fileformat=VCFv4.2\\n", encoding="utf-8")
""",
    )
    _make_exec(
        Path(tools["kanpig"]),
        """#!/usr/bin/env python3
import sys
from pathlib import Path
args = sys.argv[1:]
out = Path(args[args.index("--out") + 1])
rnames = Path(args[args.index("--rnames") + 1])
samples = []
for i, arg in enumerate(args):
    if arg == "--sample":
        samples.append(args[i + 1])
out.parent.mkdir(parents=True, exist_ok=True)
rnames.parent.mkdir(parents=True, exist_ok=True)
out.write_text(
    "##fileformat=VCFv4.2\\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\\"GT\\">\\n"
    "##FORMAT=<ID=FT,Number=1,Type=Integer,Description=\\"FT\\">\\n"
    "##FORMAT=<ID=SQ,Number=1,Type=Integer,Description=\\"SQ\\">\\n"
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\\"GQ\\">\\n"
    "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\\"PS\\">\\n"
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\\"DP\\">\\n"
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\\"AD\\">\\n"
    "##FORMAT=<ID=KS,Number=.,Type=Integer,Description=\\"KS\\">\\n"
    f"#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t{samples[0]}\\t{samples[1]}\\n"
    "chr1\\t10\\tsv1\\tA\\t<DEL>\\t60\\tPASS\\t.\\tGT:FT:SQ:GQ:PS:DP:AD:KS\\t0|1:0:50:50:1:10:0,10:99\\t0|0:0:10:10:1:12:12,0:99\\n"
    "chr1\\t20\\tsv2\\tA\\t<DEL>\\t60\\tPASS\\t.\\tGT:FT:SQ:GQ:PS:DP:AD:KS\\t0|0:0:10:10:1:12:12,0:99\\t1|0:0:50:50:1:11:0,11:99\\n"
    "chr1\\t30\\tsv3\\tA\\t<DEL>\\t60\\tPASS\\t.\\tGT:FT:SQ:GQ:PS:DP:AD:KS\\t0|1:0:50:50:1:10:4,6:99\\t0|1:0:50:50:1:10:5,5:99\\n",
    encoding="utf-8"
)
rnames.write_text("sv_id\\tread_name\\n", encoding="utf-8")
""",
    )
    return split_dir, reference, tools


class TestFilterSampleSpecificByAd(unittest.TestCase):
    def test_filter_sample_specific_by_ad(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            vcf = root / "merged.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n"
                "##FORMAT=<ID=FT,Number=1,Type=Integer,Description=\"FT\">\n"
                "##FORMAT=<ID=SQ,Number=1,Type=Integer,Description=\"SQ\">\n"
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GQ\">\n"
                "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"PS\">\n"
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n"
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">\n"
                "##FORMAT=<ID=KS,Number=.,Type=Integer,Description=\"KS\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
                "chr1\t10\tsv1\tA\t<DEL>\t60\tPASS\t.\tGT:FT:SQ:GQ:PS:DP:AD:KS\t0|1:0:0:0:1:12:6,6:99\t0|0:0:0:0:1:15:15,0:99\n"
                "chr1\t20\tsv2\tA\t<DEL>\t60\tPASS\t.\tGT:FT:SQ:GQ:PS:DP:AD:KS\t0|0:0:0:0:1:15:15,0:99\t1|0:0:0:0:1:12:3,9:99\n"
                "chr1\t30\tsv3\tA\t<DEL>\t60\tPASS\t.\tGT:FT:SQ:GQ:PS:DP:AD:KS\t0|1:0:0:0:1:12:5,7:99\t0|1:0:0:0:1:12:5,7:99\n",
                encoding="utf-8",
            )
            import subprocess
            subprocess.run(["bgzip", "-f", str(vcf)], check=True)
            subprocess.run(["tabix", "-f", "-p", "vcf", str(vcf) + ".gz"], check=True)
            outputs = _filter_sample_specific_by_ad(
                kanpig_vcf_gz=Path(str(vcf) + ".gz"),
                sample_a_label="A",
                sample_b_label="B",
                min_total_ad=5,
                min_target_alt_ad=1,
                other_max_alt_ad=0,
                output_dir=root / "out",
            )
            self.assertTrue(outputs["sample_a_only"].exists())
            self.assertTrue(outputs["sample_b_only"].exists())
            self.assertTrue(outputs["shared"].exists())


class TestSvPostProcessingMain(unittest.TestCase):
    def test_main_runs_end_to_end(self):
        with tempfile.TemporaryDirectory() as td:
            split_dir, reference, tools = _build_two_group_env(Path(td))
            summary = sv_post_processing_main(
                [
                    "--split-dir", str(split_dir),
                    "--reference", str(reference),
                    "--groups", "A,B",
                    "--bcftools-bin", tools["bcftools"],
                    "--truvari-bin", tools["truvari"],
                    "--kanpig-bin", tools["kanpig"],
                    "--output-dir", str(Path(td) / "run"),
                ]
            )
            self.assertTrue(Path(summary["collapsed_sorted_vcf"]).exists())
            self.assertTrue(Path(summary["kanpig_merged_vcf"]).exists())
            self.assertTrue(Path(summary["sample_a_only_vcf"]).exists())
            self.assertTrue(Path(summary["sample_b_only_vcf"]).exists())
            self.assertTrue(Path(summary["shared_vcf"]).exists())

    def test_resolve_args_requires_two_groups(self):
        parser = type("Args", (), {
            "split_dir": "/tmp/x",
            "reference": "/tmp/r",
            "groups": "A",
            "output_dir": None,
            "sample_id": None,
            "mosaic_filter_expression": "INFO/MOSAIC=1",
            "min_total_ad": 5,
            "min_target_alt_ad": 1,
            "other_max_alt_ad": 0,
            "bcftools_bin": "/bin/true",
            "truvari_bin": "/bin/true",
            "kanpig_bin": "/bin/true",
            "threads": 16,
            "kanpig_seqsim": 0.8,
            "kanpig_sizesim": 0.85,
        })
        with self.assertRaises(ValueError):
            _resolve_args(parser)
