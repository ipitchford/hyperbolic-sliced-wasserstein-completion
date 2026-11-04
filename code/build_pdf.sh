#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_dir"
mkdir -p output
pandoc manuscript.md --standalone --pdf-engine=xelatex --resource-path=. \
  --output output/hyperbolic_slicing_manuscript.pdf
pandoc manuscript.md --standalone --to=latex --resource-path=. \
  --output output/hyperbolic_slicing_manuscript.tex
