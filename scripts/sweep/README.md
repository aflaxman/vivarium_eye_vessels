# Calibration sweep tooling

The scripts the realism rounds (`docs/realism_roadmap.md`, passes five
onward) use to evaluate model-spec candidates across random seeds. They
wrap `vivarium_eye_vessels.vnv.calibrate`; the conventions are the ones the
roadmap's V&V section describes.

Run everything from the repository root with the package installed
(`pip install -e .`). Outputs go to `$SWEEP_DIR` (default `sweep_out/`,
git-ignored by convention; set the variable to keep sweeps elsewhere).

## Workflow of a round

1. **Job file.** One candidate per line: a short name, a tab, and a JSON
   object of dotted spec keys to override, e.g.

   ```
   vel018	{"particles.terminal_velocity": 0.18}
   v18_si18	{"particles.terminal_velocity": 0.18, "path_splitter.split_interval": 18}
   ```

   Anything not overridden comes from `model_spec.yaml` *as it is on disk
   when the run starts*, so edit the spec only when no queued run still
   depends on the old values.

2. **Sweep.** `scripts/sweep/sweep_jobs.sh JOBS.tsv [SEEDS] [PARALLEL]`
   runs every name x seed as its own process, `PARALLEL` at a time (one run
   is one core and 10-20 minutes; oversubscribing cores only slows every
   run). Default seeds are the five calibration seeds `7,42,909,2024,123456`;
   the later rounds added `31,77,5150`. A log ending in `MULTI-DONE` is
   complete and is skipped on relaunch, so an interrupted sweep can be
   restarted with the same command. Probe a knob on the spec seed first
   (`SEEDS=123456`), then take the survivors to all seeds.

3. **Read.** `python scripts/sweep/seed_matrix.py NAME ...` prints the
   per-seed score with paired / arterial / venous perfusion;
   `python scripts/sweep/sweep_table.py NAME ...` prints the mean score
   decomposition by term and the mean statistics. Pick the candidate with
   the best mean *and no seed materially worse*; when two are within
   noise, prefer the one with no new code.

4. **Gate.** `vnv_contact_sheet CANDIDATE.yaml --output-dir DIR` on the
   held-out seeds (11/202/909/4242, never used to choose knobs): every
   seed must colonize >= 95% of the tissue. Then `vnv_compare` and
   `vnv_growth_gif` regenerate `docs/vnv/` for the PR.

## Diagnostics

- `run_sim_dump.py OUT.pkl '{overrides}' [SEED] [STEPS]` pickles a
  finished run (particles, edges, superficial raster, geometry) so
  estimators can be developed against it without re-simulating.
- `score_dump.py DUMP.pkl ...` re-scores pickles with the current
  `calibrate.TARGETS` and lists the terms that cost half a point or more.
- `coverage_course.py '{overrides}' [SEED]` prints the growth time course:
  per-tree coverage, paired perfusion, live tips and where tips die.

## Seeds and honesty

Calibration seeds choose knobs; held-out seeds only confirm them. Report
negative results (knobs swept and rejected) in the CHANGELOG and roadmap
with their scores, and revert code that does not ship rather than leaving
it default-off.
