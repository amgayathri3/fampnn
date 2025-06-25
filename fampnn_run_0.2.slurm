#!/bin/bash
#SBATCH -p pbatch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH -t 1-00:00:00
#SBATCH -o logs/fampnn/fampnn%A_%a.out
#SBATCH -e logs/fampnn/fampnn%A_%a.err
#SBATCH --array=0-5

INPUT_DIRS=(
    6a3w_groove_cleaned
    6a3w_knob_cleaned
    6my4_groove_cleaned
    6my4_knob_cleaned
    8z8m_groove_cleaned
    8z8m_knob_cleaned
)

INP=${INPUT_DIRS[$SLURM_ARRAY_TASK_ID]}
# Use absolute paths
INP_PATH=/usr/workspace/ayyarmanjula1/RFantibody/Files/Cleaned_BigRun/$INP
OUTP=/usr/workspace/ayyarmanjula1/RFantibody/Files/Cleaned_BigRun/FAMPNN_Output_0.2/${INP}_fampnn

mkdir -p "$OUTP"

source ~/.bashrc
conda activate fampnn
cd /usr/workspace/ayyarmanjula1/fampnn/


python run_fampnn_wrapper.py --pdb_dir "$INP_PATH" --out_dir "$OUTP" --temperature 0.2 --num_samples 3
