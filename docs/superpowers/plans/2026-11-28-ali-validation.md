# Ali footprint and SNR validation implementation plan

> Execution: continue directly as requested, using cached products from 5c8cd56; no new injections, depth reconstruction, atmospheres or Galactic-model changes.

**Goal:** resolve footprint and SNR/domain definitions and execute the remaining count validation.
**Architecture:** retain completed maps and recovery distributions; add only footprint reconciliation and radial-domain diagnostics, with explicit scientific acceptance checks.
**Tech stack:** existing Python/NumPy/SciPy/Astropy environment and GitHub connector.
**Spec:** user's follow-up to commit 5c8cd56.

## Constraints and review focus

- Do not tune an area cut or distance cutoff to historical scalars.
- Preserve the density functions, atmosphere reaches and completed recovery products.
- Distinguish geometric coverage, usable masked coverage, source detection and photometric cuts.
- A real source can scatter across an intrinsic flux boundary; no hard boundary is implied by a measured SNR cut.
- Check the zero-source selection limit and radial convergence before certifying a count.

## Tasks

- [x] Trace catalog methodology in primary sources and reconcile empirical versus native SNR using cached measurements.
- [x] Build the footprint accounting ledger; quantify existing star masks and coverage/edge definitions on cached native masks. Exact historical mask attribution remains unidentified and is explicitly reported.
- [x] Integrate radial shells with the native SNR convention, reuse the F814W empirical probability and quantify flux-boundary scatter and the noise floor.
- [x] Test limiting cases, convergence and available uncertainties. Acceptance fails: the all-distance integral diverges as R^0.23 under the unchanged halo; no finite final count is adopted.
- [x] Save source-backed conclusions, numerical tables, figures and defensible Section 4.6 limitation wording; obtain independent scientific review.
- [ ] Commit the reviewed validation resolution to GitHub and verify the branch head.

Required scientific dependency: a validated joint source-associated catalog detection/selection kernel with a sufficiently vanishing zero-flux tail, or an independently justified finite population domain. Do not silently invent either. Completed F814W processing remains reused.
