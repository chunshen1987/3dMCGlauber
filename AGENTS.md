# Repository guidance

3dMCGlauber is a C++17 Monte Carlo Glauber code for modeling high-energy nuclear collisions. Use the CMake targets and the executable's normal entry points to verify changes.

## Scientific and numerical work

- Preserve the model's equations, units, coordinate conventions, and observable definitions unless the task calls for a physics change. State any new approximation and its range of validity. Distinguish momentum rapidity, pseudorapidity, and spacetime rapidity when they matter to the observable.
- Record the projectile and target, center-of-mass energy, event selection, centrality estimator and range, physics options, and input table versions relevant to a result. Keep distance, time, energy, and cross-section units explicit.
- Use fixed seeds to check reproducibility. For changes intended to preserve the simulation, compare the random-number stream and event-level outputs where practical. If a change intentionally alters the sampling, compare the relevant distributions across multiple seeds and report statistical uncertainties.
- Check relevant physics limits and invariants, such as probability bounds, symmetry, or conservation laws. Choose numerical tolerances appropriate to the observable and explain them when they affect the conclusion.
- Compare performance builds with equivalent workloads and compiler settings. Record the baseline revision, build flags, timing method, hardware, and attempted and accepted event counts where applicable. Verify relevant physics outputs alongside performance measurements.
- Use an out-of-source CMake build and run simulations in a separate directory so generated data stays out of the source tree. Use existing tests and focused simulation checks for the behavior under study.

## Investigation reports

- Save a Markdown report under `doc/` for every future investigation, including investigations that produce only findings or end inconclusively. Include the report in Git with the related changes, or in a documentation commit when the report is the only deliverable.
- Record the question, baseline revision and local changes, reproducible commands and inputs, results, validation, and limitations. Cite relevant source files or primary physics references. Keep the report usable after temporary build and output directories are gone.
- Use the filename and outline in [doc/README.md](doc/README.md).

## Git commits

- Give each commit a concise subject and a brief body summarizing what changed and why. Include the most relevant verification when available.
- End each AI-authored commit message with `Signed-off-by: Codex (AI agent) <codex@localhost>`. If another AI agent makes the commit, use that agent's name and sign-off identity. Keep the configured Git author and committer identities.
