set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
Rscript render_quarto_report.R qc test_datasets/Human_Kidney/parameters/human_kidney_test.yaml

Rscript render_quarto_report.R de test_datasets/Human_Kidney/parameters/human_kidney_test.yaml
