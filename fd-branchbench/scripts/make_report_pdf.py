#!/usr/bin/env python
"""Render report/index.html (the published v3 report) to a print PDF.

Usage: make_report_pdf.py [report/index.html] [../folddisco_3.0_report.pdf]
Needs weasyprint in the venv: .venv/bin/pip install weasyprint
"""
import pathlib, sys
from weasyprint import HTML

HERE = pathlib.Path(__file__).resolve().parents[1]
SRC = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "report/index.html"
OUT = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else HERE.parent / "folddisco_3.0_report.pdf"

PRINT = """
@page {
  size: A4; margin: 17mm 15mm 16mm;
  @bottom-center { content: counter(page) " / " counter(pages);
    font-family: "IBM Plex Mono", monospace; font-size: 8pt; color: #7c8186; }
  @bottom-left { content: "Folddisco 3.0 - update report";
    font-family: "IBM Plex Sans", sans-serif; font-size: 8pt; color: #7c8186; }
}
body { background: #fff; font-size: 10pt; line-height: 1.5; }
.wrap { max-width: none; padding-inline: 0; padding-block: 0; }
.measure { max-width: none; }
header.top { padding-block: 0 18px; margin-bottom: 20px; }
h1 { font-size: 24pt; }
.lede { font-size: 11.5pt; }
h2 { font-size: 15pt; margin: 22px 0 4px; padding-top: 12px; break-after: avoid; break-before: auto; }
h3 { font-size: 11.5pt; margin: 16px 0 3px; break-after: avoid; }
p, li { orphans: 2; widows: 2; }
figure, .tablewrap, pre, .tiles, .note { break-inside: avoid; }
figure { margin: 14px 0; }
figure img { max-height: 105mm; width: auto; max-width: 100%; margin: 0 auto; }
.tablewrap { overflow: visible; }
table { font-size: 8pt; table-layout: auto; }
th, td { padding: 5px 8px; white-space: normal; }
td.n, th.n { white-space: nowrap; }
.swatches { gap: 4px 14px; }
pre { font-size: 7.6pt; line-height: 1.4; white-space: pre-wrap; }
code { font-size: 0.86em; }
.tile .v { font-size: 15pt; }
footer { margin-top: 24px; }
a { color: #c2286b; text-decoration: none; }
"""

html = SRC.read_text()
head, rest = html.split("<style>", 1)
css, body = rest.split("</style>", 1)
# appended to the page's own stylesheet: a stylesheet passed to write_pdf is user-level and
# would lose to these author rules.
doc = (f'<!doctype html><html lang="en"><head><meta charset="utf-8">{head}'
       f'<style>{css}\n{PRINT}</style></head><body>{body}</body></html>')
HTML(string=doc, base_url=str(SRC.parent.resolve())).write_pdf(OUT)
print("wrote", OUT, f"{OUT.stat().st_size/1e6:.1f} MB")
