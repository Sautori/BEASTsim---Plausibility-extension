# Extending BEASTsim

Code for my MSc thesis, *Extending BEASTsim with multi-task benchmarking and biological
plausibility assessment*. The project adds a set of task metrics and a candidate biological
plausibility layer on top of BEASTsim. It runs them on five real spatial transcriptomics
datasets across three platforms, then reads the output through several plausibility criteria
instead of one score.

This repository holds two things: the evaluation pipeline that produced the run, and the
plotting code that turns the run outputs into every data figure in the thesis.

## Layout

```
pipeline/   the evaluation: metric streams, copy-risk, the plausibility module, the manifest
figures/    one script per figure, plus make_figures.py to build them all
ucloud/     the batch script for the full run on SDU UCloud
data/       the run summary tables the figures read. h5ad tissues come from Zenodo
```

## Data

The processed datasets and the simulated outputs are on Zenodo:
https://zenodo.org/records/15775564

The summary tables the figures read are already in `data/`:

- `metric_long.csv` - every metric value for every (dataset, method-variant) pair
- `component_evidence.csv` - the four plausibility components per pair
- `named_plausibility_baseline.csv` - the named candidate summaries and copy-risk
- `expected_vs_observed.csv` - the expected direction of each perturbation control
- `perturbation_response.csv` - the measured response of each control
- `spatialsimbench_input_manifest.csv` - which source file each pair was read from

The per-pair `.h5ad` tissues are too large to keep here. They are in the Zenodo deposit.
Download them into `data/` if you want the tissue figures (the gallery, the real-versus-simulated
panel, the neighbour matrices, the per-cell-type overlap). The manifest lists the filenames the
scripts expect.

## Run on your own computer

The figures build from the summary tables already in `data/`, so this part needs no downloads:

```
pip install -r requirements.txt
cd figures
python make_figures.py
```

Figures land in `outputs/`. The tissue figures are skipped unless the `.h5ad` files from Zenodo
are sitting in `data/`.

To redo the whole evaluation locally, not just the figures, put the Zenodo `.h5ad` datasets in
`data/` and run:

```
cd pipeline
python thesis_pipeline.py --manifest ../data/spatialsimbench_input_manifest.csv --data-root .. --output-root ../outputs
```

That writes `metric_long.csv` and the other tables the figures read. Each pair peaks at about
1.4 to 3.3 GB, so keep `--workers` low (2 to 4) to stop the parallel pairs stacking past your
memory. Lower `--corr-sample` and the cell caps to trade a little accuracy for time. The manifest
is built by `pipeline/build_beastsim_faithful_manifest.py`, the metric definitions live in
`pipeline/spatialsimbench_metrics.py`, and the plausibility components and named summaries in
`pipeline/bp_module.py`.

### Replicating the main comparison

By default the run uses five seeds (`SEEDS="0,1,2,3,4"`) so it carries run-to-run uncertainty on
the main comparison. Set `SEEDS="0"` for the original single pass. The pipeline evaluates every
pair under each seed, tags each output row with its `seed`, and `aggregate_seeds.py` writes the
replicate statistics: per-variant Balanced Realism mean and spread (`seed_stability_per_variant.csv`)
and a sign test of scDesign3 over reference-free SRTsim on the three shared datasets
(`seed_stability_shared_comparison.csv`). Seeds vary the evaluation stochasticity (clustering,
subsampling) on the fixed simulated outputs, so this measures whether a margin survives evaluation
noise. It does not re-draw the simulators, which would need their generators rather than the
deposited outputs.

## Run on UCloud

The full run was one batch job on SDU UCloud. Everything ships in a single self-contained zip
that the batch script unwraps for you:

1. From the repo root, build the bundle: `ucloud/make_upload.sh <folder-of-h5ad-files>`. It writes
   one `ucloud_bundle.zip` holding `scripts/`, `manifests/`, `data/`, `requirements.txt`, and
   `run_all.sh`. The `.h5ad` datasets come from the Zenodo deposit. Leave the data argument out and
   it leaves a note in `data/` so you can add them before zipping.
2. Upload that one `ucloud_bundle.zip` to a UCloud drive. Nothing else.
3. Start the **terminal-ubuntu** app, mount the drive, and attach `run_all.sh` as the batch script.
   `run_all.sh` finds the zip, unwraps it, installs the requirements, and runs.
4. Outputs land in `outputs/` on the drive: `metric_long.csv`, the component and summary tables,
   the replicate statistics, an environment record, and the run log.

Workers default to the core count. Set `WORKERS` to override it, and lower `--corr-sample` and
`--pert-max-obs` to run faster on a small machine.

## Compute

The run used **cpu-amd-zen5**: one node, 128 cores, 384 GB of RAM, a two hour limit. That size is
not required. The work is parallel across the 29 pairs and finishes on far smaller hardware. The
node reflects what was on hand, not what the job needs.

## Citation

If you use this code, please cite the thesis and the Zenodo deposit.

## License

MIT, see `LICENSE`.
