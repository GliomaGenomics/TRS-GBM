#!/bin/bash
#SBATCH --job-name=trs_single_cell
#SBATCH --time=05:00:00
#SBATCH --mem=240G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=hpc/logs/%x_%j.out

# Every analysis that reads one single-cell dataset, run on the real cells: the
# pseudobulks first, then each panel drawn from the cells or from those
# pseudobulks.
#
# The job runs in its own folder on scratch, OUTPUT_ROOT/<dataset>_<job id>,
# holding a copy of the code, the config, the pixi files and data/ (whose
# single-cell files are links). The scripts write where they always do relative
# to that copy, so the figures, source data and pseudobulks land there, beside
# the code and config that made them, and nothing large is written to the home
# directory. The pixi environment cannot move, so the copied scripts run in the
# environment installed in the clone.
#
# Once, in the clone in the home directory, on a login node:
#     pixi install
#     mkdir -p hpc/logs
#     ln -s <real h5ad on scratch> data/single-cell/<file named in config/datasets.yaml>
#
# Submit from the root of that clone, one job per dataset:
#     sbatch hpc/run_single_cell.sh nomura_2025
#     sbatch hpc/run_single_cell.sh wang_2022

set -euo pipefail

OUTPUT_ROOT=$SCRATCH/TRS

dataset=$1
repository=$PWD
run=$OUTPUT_ROOT/${dataset}_${SLURM_JOB_ID}

in_environment() {
    pixi run --manifest-path "$repository/pixi.toml" "$@"
}

# Outputs drawn from the 1/100 dummy subsample would pass for the real ones.
cells=$(in_environment python -c "
import sys
sys.path.insert(0, 'code')
from common import cohort
print(cohort('single_cell', sys.argv[1])['path'])
" "$dataset")
if [[ $(readlink -f "$cells") == *dummy* ]]; then
    echo "$cells is the dummy file $(readlink -f "$cells")" >&2
    exit 1
fi

mkdir -p "$run"
cp -a code config data pixi.toml pixi.lock "$run/"
cd "$run"

in_environment python code/pseudobulk.py "$dataset"
in_environment python code/composition_vs_axis_score.py "$dataset"
in_environment python code/axis_score_composition_regression.py "$dataset" --component celltype
in_environment python code/axis_score_composition_regression.py "$dataset" --component state
in_environment python code/celltype_axis_score_distribution.py "$dataset"
in_environment python code/malignant_state_axis_score_by_stage.py "$dataset"
in_environment python code/single_cell_patient_table.py "$dataset"
