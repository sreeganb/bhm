#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd
#$ -l scratch=200G
#$ -l h_rt=24:00:00

# ————————————————————
# 0. Project & scratch
# ————————————————————
PROJECT_ROOT="/wynton/home/sali/sreeganb/git/toy_model"

# Use the scheduler’s TMPDIR if set; otherwise fall back.
if [[ -z "$TMPDIR" ]]; then
  echo "WARNING: \$TMPDIR not set by SGE; falling back to /tmp/\$USER"
  TMPDIR="/tmp/$USER/job_$JOB_ID"
  mkdir -p "$TMPDIR"
fi
echo "Using scratch dir: $TMPDIR"
cd "$TMPDIR" || exit 1

# ————————————————————
# 1. Copy only what you need
# ————————————————————
# Fix PROJECT_ROOT if it was wrong
if [[ ! -d "$PROJECT_ROOT" ]]; then
  echo "ERROR: PROJECT_ROOT not found: $PROJECT_ROOT" >&2
  exit 1
fi

# Create the local subdirs
mkdir -p output_analysis results

# Use rsync (fast, incremental, preserves perms) to pull down
rsync -av --delete "$PROJECT_ROOT/output_analysis/" "$TMPDIR/output_analysis/"
rsync -av --delete "$PROJECT_ROOT/results/"         "$TMPDIR/results/"

# Copy scripts
cp "$PROJECT_ROOT/"*.py "$TMPDIR/"

# ————————————————————
# 2. Run the sampler
# ————————————————————
echo "Starting MCMC at $(date)"
python3.11 run_samplers.py 2> output_tetramer.out
echo "Finished MCMC at $(date)"

# ————————————————————
# 3. Push back outputs
# ————————————————————
# Backup any old output_analysis
if [[ -d "$PROJECT_ROOT/results/output_analysis" ]]; then
  mv "$PROJECT_ROOT/results/output_analysis" \
     "$PROJECT_ROOT/results/output_analysis.$(date +%Y%m%d_%H%M%S)"
fi

# Mirror results back
rsync -av --delete "$TMPDIR/output_analysis/" "$PROJECT_ROOT/results/output_analysis/"
cp output_tetramer.out "$PROJECT_ROOT/results/"

# ————————————————————
# 4. Job summary
# ————————————————————
echo "Job complete."
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
