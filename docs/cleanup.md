# Revision and preservation record

This revision reorganizes the core physical hierarchy into M0 ideal inverse square, M1 finite geometry with a shared offset, and M2 geometry plus effective attenuation. Every model now shares the same background/dead-time response and statistically comparable nuisance treatment. The new M0 was refitted; it is not the old uncorrected M0.

Before changing files, the complete tracked previous src/, results/, docs/, README, requirements and reproduction entry point were copied to archive/2026-11-03_M0-M3/. All previous numerical results and figures are retained there, explicitly labelled historical. Root results/ now holds only current generated outputs. Historical report source and images remain in report/ unchanged. Raw CSV bytes and provenance hashes are unchanged.

The existing validated geometry integrator was reused without changing its algorithm or constants; unused imports were removed. The optimizer now applies the same physical-positive latent-distance domain to all models rather than inconsistent per-model sigma clipping. Broad remaining numerical bounds are checked inactive. This makes the comparison consistent and changes no apparatus inputs.

P1 preserves the previous shifted point-source/small-detector geometry; P2 adds attenuation to test its effect on μ. These are sensitivity checks, not additional core model stages. The final three-panel plotting routine is separate from reusable models and comparison code. No redundant historical output was deleted; archived source and results are independent of root output regeneration.

Current numerical/figure verification is documented in docs/verification.md. Original consolidation details remain in archive/2026-11-03_M0-M3/docs/cleanup.md. The original project under 2026-10-29/r was not written to in either revision workflow.
