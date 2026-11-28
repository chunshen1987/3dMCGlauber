# Minimum-bias centrality table reuse

## Question and scope

Can later runs reuse the 50,000-event minimum-bias sample used to set
centrality cuts on the number of produced strings, `Nstrings`? The table stores
the descending list of string counts. The calibration sampling and percentile
indexing are in [EventGenerator.cpp](../src/EventGenerator.cpp). The production
trigger applies inclusive `Nstrings` and participant-count limits. The change
introduces no physics approximation.

## Setup and reproduction

The baseline was `6cf5ea61ca5ffe2e51bff6365da18e179a9454be` with a clean worktree.
The comparison used Linux x86_64, GNU C++ 15.2.0, and the CMake defaults
`-g -O3 -Wall -std=c++17`. The updated source also used `-Dunittest=ON`
to build the test executables.

From the updated repository, build both revisions and prepare isolated run
directories with the same checked input and quark tables.

```sh
repo=$PWD
scratch=$(mktemp -d)
mkdir -p "$scratch/baseline" "$scratch/base-run" "$scratch/new-run"
git archive 6cf5ea61ca5ffe2e51bff6365da18e179a9454be | tar -x -C "$scratch/baseline"
cmake -S "$scratch/baseline" -B "$scratch/base-build"
cmake -S "$repo" -B "$scratch/new-build" -Dunittest=ON
cmake --build "$scratch/base-build" --target 3dMCGlb.e -j 8
cmake --build "$scratch/new-build" -j 8
ln -s "$repo/tables" "$scratch/base-run/"
ln -s "$repo/tables" "$scratch/new-run/"

run_oxygen() {
    "$1" 20 "$repo/input" 12345 Projectile=O Target=O \
        only_event_statistics=1 cenMin=0 cenMax=20
}

cd "$scratch/base-run"
run_oxygen "$scratch/base-build/src/3dMCGlb.e" > baseline.log
cd "$scratch/new-run"
run_oxygen "$scratch/new-build/src/3dMCGlb.e" > cold.log
cmp events_summary.dat ../base-run/events_summary.dat
run_oxygen "$scratch/new-build/src/3dMCGlb.e" > warm.log
cp events_summary.dat warm-events.dat
run_oxygen "$scratch/new-build/src/3dMCGlb.e" > repeat.log
cmp events_summary.dat warm-events.dat
grep 'centrality cut' ../base-run/baseline.log cold.log warm.log repeat.log
wc -l minimum_bias_collisions.dat
```

This selects 20 O+O events at 200 GeV with 0–20% `Nstrings` centrality
and `2 <= Npart <= 500`. The inelastic NN cross section is 40.8344 mb.
The production impact-parameter range is 0–20 fm. Calibration accepts events
with at least two participants and starts at 0–25 fm before applying the
existing adaptive upper limit. The remaining [input](../input) settings
include sampled nuclear configurations, a minimum nucleon spacing of 0.9 fm,
valence quarks at `Q2=1` GeV², string-number fluctuations, and
`shadowing_factor=0.6`.

The quark sample files used in the O+O and p+p runs had these SHA-256 digests.

| Input file under `tables/` | SHA-256 |
| --- | --- |
| `proton_valence_quark_samples.dat` | `209996edb93ece23937fdaba90ce8040a549ac90ad01b9123448d61e99f40725` |
| `neutron_valence_quark_samples.dat` | `a688a88a98774540d1fb6f6a28f31eff004ea480013002ba2b067d0c9e17f66d` |

## Results

The baseline and updated O+O runs gave the same inclusive cuts
`16 <= Nstrings <= 41`. The baseline and updated first-run event summaries
matched byte for byte. The two runs using the saved table also matched byte
for byte. Loading preserved the table contents and modification time and
skipped the 50,000-event calibration.

Additional p+p runs used the same options and input files with
`Projectile=p Target=p minBiasTableFile=/path/to/custom.dat`. A run in a
second directory loaded the first run's table with seed 42 and 20–50%
centrality. Truncating the file to its first ten lines caused regeneration,
recreating the original bytes at seed 12345. Changing `roots` from 200 to
201 GeV replaced the incompatible table, and repeating at 201 GeV loaded
it. Pointing `minBiasTableFile` into a missing directory produced a write
warning and the requested events. A 0–100% run skipped table generation.

## Validation and limitations

The saved O+O table contained 50,000 positive, descending integer counts.
Its entries at zero-based indices 0 and 10,000 matched the reported 0% and
20% cut values. The 20 events in each O+O summary satisfied the participant,
string-count, and impact-parameter cuts, and proton plus neutron participants
summed to `Npart`. Seven of nine existing unit-test executables passed.
`unittest_Nucleus.e` segfaulted and `unittest_Quark.e` aborted at a
mass–momentum assertion on both revisions with the same test fixtures.

The compatibility check records the effective collision energy and the
resolved input settings apart from the run seed, production cuts, output
options, and table path. The external nuclear and quark table contents are
outside the check. Use a fresh calibration after changing those files.
The first run's production stream starts after calibration draws, while
reusing a table skips those draws. Report the table and the production seed
when reproducing a run. This check compared cases at fixed seeds and used
no timing benchmark.
