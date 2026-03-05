import json
import re
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd

from sniffcell.report.igvreport import (
    _customize_igvreport_html,
    _decode_igv_data_uri_json,
    _encode_igv_data_uri_json,
    build_igvreport_cli_command,
    render_igvreport_bundle,
)


class TestIgvReportHelpers(unittest.TestCase):
    def test_customize_igvreport_html_sets_basemod_and_support_highlights(self):
        with tempfile.TemporaryDirectory() as td:
            html_path = Path(td) / "index.html"
            session_uri = _encode_igv_data_uri_json(
                {
                    "locus": "chr1:100-200",
                    "reference": {"fastaURL": "data:application/gzip;base64,ZmFrZQ=="},
                    "tracks": [
                        {"type": "alignment", "name": "reads", "url": "data:application/gzip;base64,ZmFrZQ=="},
                        {"type": "variant", "name": "sv", "url": "data:application/gzip;base64,ZmFrZQ=="},
                    ],
                }
            )
            html_path.write_text(
                (
                    "<html><body><div id=\"igvDiv\"></div><script type=\"text/javascript\">\n"
                    "const tableJson = {\"headers\":[\"unique_id\"],\"rows\":[[0]]}\n"
                    f"const sessionDictionary = {{\"0\": {json.dumps(session_uri)}}}\n"
                    "let igvBrowser\n"
                    "document.addEventListener(\"DOMContentLoaded\", function () {\n"
                    "    initIGV()\n"
                    "})\n"
                    "function initIGV() {\n"
                    "    const igvDiv = document.getElementById(\"igvDiv\")\n"
                    "    const options = {\n"
                    "        sessionURL: sessionDictionary[\"0\"],\n"
                    "        showChromosomeWidget: false,\n"
                    "        showCenterGuide: true,\n"
                    "        search: false\n"
                    "    }\n"
                    "    igv.createBrowser(igvDiv, options)\n"
                    "        .then(function (b) {\n"
                    "                igvBrowser = b\n"
                    "                initTable()\n"
                    "        })\n"
                    "}\n"
                    "function initTable() {\n"
                    "    const uniqueId = \"0\"\n"
                    "    const session = sessionDictionary[uniqueId]\n"
                    "                igvBrowser.loadSession({\n"
                    "                    url: session\n"
                    "                })\n"
                    "}\n"
                    "</script></body></html>\n"
                ),
                encoding="utf-8",
            )

            _customize_igvreport_html(html_path, supporting_reads_by_session={"0": ["readA", "readB"]})

            html_text = html_path.read_text(encoding="utf-8")
            self.assertIn("sniffcellSupportingReadsBySession", html_text)
            self.assertIn("\"readA\", \"readB\"", html_text)
            self.assertIn("setHighlightedReads", html_text)
            self.assertIn("sniffcellApplySupportingReadHighlights(\"0\")", html_text)
            self.assertIn("sniffcellEnsureSupportToolbar", html_text)
            self.assertIn("Highlight supporting reads", html_text)
            self.assertIn("true read selection/filtering is not exposed by IGV.js here", html_text)

            match = re.search(r"const sessionDictionary = (.+?)\n\s*let igvBrowser", html_text, flags=re.DOTALL)
            self.assertIsNotNone(match)
            session_dict = json.loads(match.group(1))
            session_payload = _decode_igv_data_uri_json(session_dict["0"])
            self.assertEqual(session_payload["tracks"][0]["colorBy"], "basemod2")
            self.assertNotIn("colorBy", session_payload["tracks"][1])

    def test_build_igvreport_cli_command_uses_generic_site_columns(self):
        cmd = build_igvreport_cli_command(
            runner_prefix=["create_report"],
            sites_path=Path("/tmp/sites.tsv"),
            output_html=Path("/tmp/index.html"),
            reference_path="/tmp/ref.fa",
            tracks=["/tmp/a.bam", "/tmp/in.vcf.gz"],
            info_columns=["sv_locus", "primary_celltype"],
            flanking=5000,
            header_path=Path("/tmp/header.html"),
        )
        self.assertEqual(cmd[0], "create_report")
        self.assertIn("--output", cmd)
        self.assertIn("--sequence", cmd)
        self.assertIn("1", cmd)
        self.assertIn("--begin", cmd)
        self.assertIn("2", cmd)
        self.assertIn("--end", cmd)
        self.assertIn("3", cmd)
        self.assertIn("--zero_based", cmd)
        self.assertIn("--tracks", cmd)
        self.assertIn("/tmp/a.bam", cmd)
        self.assertIn("/tmp/in.vcf.gz", cmd)
        self.assertIn("--info-columns", cmd)

    def test_render_igvreport_bundle_writes_sites_manifest_and_html(self):
        with tempfile.TemporaryDirectory() as td:
            anno_dir = Path(td) / "anno_out"
            anno_dir.mkdir(parents=True, exist_ok=True)
            (anno_dir / "anno_run_manifest.json").write_text(
                json.dumps(
                    {
                        "inputs": {
                            "bam": "/tmp/default.bam",
                            "vcf": "/tmp/input.vcf.gz",
                            "reference": "/tmp/ref.fa",
                            "bed": "/tmp/dmrs.tsv",
                        }
                    }
                ),
                encoding="utf-8",
            )
            selected_df = pd.DataFrame(
                [
                    {
                        "id": "sv1",
                        "sv_chr": "chr1",
                        "sv_pos": 101,
                        "sv_len": 250,
                        "supporting_reads": "readA,readB",
                        "primary_celltype": "A",
                        "linked_celltypes": "A|B",
                        "assigned_code": "10",
                        "majority_pct": 0.9,
                        "overlap_pct": 0.8,
                        "n_supporting": 8,
                        "n_overlapped": 6,
                    }
                ]
            )
            output_dir = anno_dir / "report" / "igvreport"

            def _fake_run(cmd, check, capture_output, text):
                html_path = Path(cmd[cmd.index("--output") + 1])
                html_path.parent.mkdir(parents=True, exist_ok=True)
                session_uri = _encode_igv_data_uri_json(
                    {
                        "locus": "chr1:100-350",
                        "reference": {"fastaURL": "data:application/gzip;base64,ZmFrZQ=="},
                        "tracks": [{"type": "alignment", "name": "reads", "url": "data:application/gzip;base64,ZmFrZQ=="}],
                    }
                )
                html_path.write_text(
                    (
                        "<html><body><div id=\"igvDiv\"></div><script type=\"text/javascript\">\n"
                        "const tableJson = {\"headers\":[\"unique_id\"],\"rows\":[[0]]}\n"
                        f"const sessionDictionary = {{\"0\": {json.dumps(session_uri)}}}\n"
                        "let igvBrowser\n"
                        "document.addEventListener(\"DOMContentLoaded\", function () {\n"
                        "    initIGV()\n"
                        "})\n"
                        "function initIGV() {\n"
                        "    const igvDiv = document.getElementById(\"igvDiv\")\n"
                        "    const options = {\n"
                        "        sessionURL: sessionDictionary[\"0\"],\n"
                        "        showChromosomeWidget: false,\n"
                        "        showCenterGuide: true,\n"
                        "        search: false\n"
                        "    }\n"
                        "    igv.createBrowser(igvDiv, options)\n"
                        "        .then(function (b) {\n"
                        "                igvBrowser = b\n"
                        "                initTable()\n"
                        "        })\n"
                        "}\n"
                        "function initTable() {\n"
                        "    const uniqueId = \"0\"\n"
                        "    const session = sessionDictionary[uniqueId]\n"
                        "                igvBrowser.loadSession({\n"
                        "                    url: session\n"
                        "                })\n"
                        "}\n"
                        "</script></body></html>\n"
                    ),
                    encoding="utf-8",
                )
                return SimpleNamespace(returncode=0, stdout="ok", stderr="")

            with patch("sniffcell.report.igvreport._resolve_igvreport_runner", return_value=["create_report"]), patch(
                "sniffcell.report.igvreport.subprocess.run",
                side_effect=_fake_run,
            ) as mock_run:
                result = render_igvreport_bundle(
                    anno_output=anno_dir,
                    selected_df=selected_df,
                    output_dir=output_dir,
                    native_report_html=anno_dir / "report" / "index.html",
                    igv_bams=["/tmp/a.bam", "/tmp/b.bam"],
                    window=6000,
                )

            self.assertEqual(result["status"], "rendered")
            self.assertIn("create_report", str(result["command"]))
            self.assertEqual(mock_run.call_count, 1)

            sites_tsv = output_dir / "selected_sv_sites.tsv"
            self.assertTrue(sites_tsv.exists())
            sites_df = pd.read_csv(sites_tsv, sep="\t")
            self.assertEqual(sites_df["chrom"].tolist(), ["chr1"])
            self.assertEqual(sites_df["start"].tolist(), [100])
            self.assertEqual(sites_df["end"].tolist(), [350])
            self.assertEqual(sites_df["sv_locus"].tolist(), ["chr1:101-350"])

            manifest_payload = json.loads((output_dir / "igvreport_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest_payload["status"], "rendered")
            self.assertEqual(manifest_payload["flanking"], 6000)
            self.assertEqual(
                manifest_payload["tracks"],
                ["/tmp/a.bam", "/tmp/b.bam", "/tmp/input.vcf.gz"],
            )
            self.assertIn("sv_locus", manifest_payload["info_columns"])
            self.assertEqual(manifest_payload["session_customization"]["alignment_color_by"], "basemod2")
            self.assertEqual(manifest_payload["session_customization"]["status"], "applied")
            self.assertEqual(manifest_payload["session_customization"]["supporting_reads_available_sessions"], 1)

            header_text = (output_dir / "sniffcell_igvreport_header.html").read_text(encoding="utf-8")
            self.assertIn("Open native SniffCell report", header_text)
            self.assertIn("DNA methylation two-color mode", header_text)
            html_text = (output_dir / "index.html").read_text(encoding="utf-8")
            self.assertIn("sniffcellSupportingReadsBySession", html_text)
            self.assertIn("setHighlightedReads", html_text)
            self.assertIn("Highlight supporting reads", html_text)
            self.assertTrue((output_dir / "index.html").exists())

    def test_render_igvreport_bundle_rerenders_existing_html(self):
        with tempfile.TemporaryDirectory() as td:
            anno_dir = Path(td) / "anno_out"
            anno_dir.mkdir(parents=True, exist_ok=True)
            (anno_dir / "anno_run_manifest.json").write_text(
                json.dumps({"inputs": {"reference": "/tmp/ref.fa", "bam": "/tmp/default.bam", "vcf": "/tmp/input.vcf.gz"}}),
                encoding="utf-8",
            )
            selected_df = pd.DataFrame([{"id": "sv1", "sv_chr": "chr1", "sv_pos": 10, "sv_len": 10}])
            output_dir = anno_dir / "report" / "igvreport"
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "index.html").write_text("<html>existing</html>", encoding="utf-8")
            (output_dir / "igvreport_manifest.json").write_text(
                json.dumps({"status": "rendered", "igvreport_command": "create_report cached"}),
                encoding="utf-8",
            )

            def _fake_run(cmd, check, capture_output, text):
                html_path = Path(cmd[cmd.index("--output") + 1])
                session_uri = _encode_igv_data_uri_json(
                    {
                        "locus": "chr1:10-20",
                        "reference": {"fastaURL": "data:application/gzip;base64,ZmFrZQ=="},
                        "tracks": [{"type": "alignment", "name": "reads", "url": "data:application/gzip;base64,ZmFrZQ=="}],
                    }
                )
                html_path.write_text(
                    (
                        "<html><body><div id=\"igvDiv\"></div><script type=\"text/javascript\">\n"
                        "const tableJson = {\"headers\":[\"unique_id\"],\"rows\":[[0]]}\n"
                        f"const sessionDictionary = {{\"0\": {json.dumps(session_uri)}}}\n"
                        "let igvBrowser\n"
                        "document.addEventListener(\"DOMContentLoaded\", function () {\n"
                        "    initIGV()\n"
                        "})\n"
                        "function initIGV() {\n"
                        "    const igvDiv = document.getElementById(\"igvDiv\")\n"
                        "    const options = {\n"
                        "        sessionURL: sessionDictionary[\"0\"],\n"
                        "        showChromosomeWidget: false,\n"
                        "        showCenterGuide: true,\n"
                        "        search: false\n"
                        "    }\n"
                        "    igv.createBrowser(igvDiv, options)\n"
                        "        .then(function (b) {\n"
                        "                igvBrowser = b\n"
                        "                initTable()\n"
                        "        })\n"
                        "}\n"
                        "function initTable() {\n"
                        "    const uniqueId = \"0\"\n"
                        "    const session = sessionDictionary[uniqueId]\n"
                        "                igvBrowser.loadSession({\n"
                        "                    url: session\n"
                        "                })\n"
                        "}\n"
                        "</script></body></html>\n"
                    ),
                    encoding="utf-8",
                )
                return SimpleNamespace(returncode=0, stdout="ok", stderr="")

            with patch("sniffcell.report.igvreport._resolve_igvreport_runner", return_value=["create_report"]), patch(
                "sniffcell.report.igvreport.subprocess.run",
                side_effect=_fake_run,
            ) as mock_run:
                result = render_igvreport_bundle(
                    anno_output=anno_dir,
                    selected_df=selected_df,
                    output_dir=output_dir,
                    native_report_html=anno_dir / "report" / "index.html",
                    igv_bams=None,
                    window=5000,
                )

            self.assertEqual(result["status"], "rendered")
            self.assertIn("create_report", str(result["command"]))
            self.assertEqual(mock_run.call_count, 1)
            html_text = (output_dir / "index.html").read_text(encoding="utf-8")
            self.assertIn("sniffcellEnsureSupportToolbar", html_text)


if __name__ == "__main__":
    unittest.main()
