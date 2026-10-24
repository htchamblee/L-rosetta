# Mixed D-RNA / L-RNA modeling (FARFAR / rna_denovo)

This note covers what this fork changes so that D-RNA and L-RNA (mirror-image RNA) can be modeled
in the same pose, and what is still D-only.

## Short version

* Write L-RNA residues as `1`, `2`, `3`, `4` (L-a, L-c, L-g, L-u) in `-sequence` or the FASTA file.
  `a c g u` are still D-RNA. The digits are input shorthand: `1` is expanded to `a[LRAD]`, the
  annotated form Rosetta already supports. Inside Rosetta the L types keep the one-letter codes
  `a c g u`, so `pose.sequence()` prints letters and the annotated sequence shows `[LRAD]` etc.
  Example for a D strand and an L strand:

  ```
  >D chain A:1-8
  ggcgaagc
  >L chain B:1-8
  33232213
  ```

* No custom params files are needed. Rosetta generates the L-RNA residue types (`LRAD`, `LRCY`,
  `LRGU`, `LURA`, PDB names ` 0A`, ` 0C`, ` 0G`, ` 0U`) from the D types with the `L_RNA` patch
  (`patches/nucleic/rna/L_RNA.txt`). Their geometry is an exact mirror of the D types.
* Everything below turns on automatically when the pose contains L-RNA. In a pure D-RNA pose
  every changed code path reduces to the original arithmetic, so results are unchanged.
* Override: `-rna:denovo:invert_torsions_for_l_rna false` disables fragment and jump mirroring
  (the old behavior). Only use it if `-vall_torsions` points at a library that is already mirrored
  and the pose is pure L-RNA. Scoring stays chirality-aware either way.
* The first time an L-RNA window receives a fragment, the `protocols.rna.denovo.fragments.FullAtomRNA_Fragments`
  tracer prints one line:
  `L-RNA detected: inserting mirror-image (sign-inverted) vall fragments into L-RNA residues ...`

## How mirroring works

Reflecting a structure keeps every bond length and angle and negates every dihedral. The vall,
the jump library, the ideal-geometry PDBs and every knowledge-based potential in Rosetta come from
D-RNA. So an L-RNA residue is handled as the mirror image of a D-RNA residue:

| Quantity | L-RNA handling |
|---|---|
| Backbone torsions α–ζ, χ, ν2, ν1, 2'-OH | negated (`rna_chirality_factor()` = −1) |
| Non-main-chain sugar coordinates (C3'/C4'/C5' stub frame) | local `z` negated (`mirror_non_main_chain_sugar_coords()`) |
| Base-pair jumps (`RNA_JumpMover`) | `M → F·M·F`, `t → F·t`, `F = diag(1,1,−1)` |
| Base-pair-step chunks (`-bps_moves`) | library coordinates reflected (x → −x) for all-L steps; mixed D/L steps skipped |
| Base coordinate frame used for scoring | `z` axis negated (`get_rna_base_coordinate_system_chirality_normalized()`) |
| Ideal-geometry reference residues (`RNA_IdealCoord`) | mirrored copies (x → −x) of the D PDBs |

The stub-based rules follow from how Rosetta builds frames: the third axis is always a cross
product, so reflecting a residue keeps the local x and y coordinates and negates z.

## What changed, by file

### Detection and helpers
* `core/chemical/rna/util.{hh,cc}`: `rna_chirality_factor()`,
  `get_rna_base_coordinate_system_chirality_normalized()` and its `_stub_` variant. `get_residue_pucker_state()`
  and `get_residue_base_state()` (syn/anti) now classify the mirrored torsion.
* Chirality comes from `ResidueType::is_l_rna()` (the `L_RNA` property).

### Fragment insertion (vall)
* `core/fragment/rna/FullAtomRNA_Fragments.{hh,cc}`
  * `insert_fragment()` negates every torsion and mirrors the sugar coordinates when the window is
    L-RNA. All callers get this: `apply_random_fragment`, `apply_best_fragment`, and the symmetry copies.
  * `window_has_uniform_rna_chirality()`: a window that mixes D and L residues is never filled.
    `insert_fragment()` skips it with a one-time warning.
  * Vall lookup maps L residues to the D sequence (`1` → `a`, ...), so homology exclusion and
    sequence matching work unchanged.
* `protocols/rna/denovo/movers/RNA_FragmentMover.cc`: mixed-chirality windows are removed from the
  insertion map, so they are never picked and no Monte Carlo moves are wasted on them.

### Base-pair jumps
* `core/import_pose/RNA_JumpMover.cc`: a jump drawn from the D-RNA jump library between two L-RNA
  residues is mirrored. This covers the initial pairings and later `random_jump_change` moves.

### Base-pair steps
* `core/import_pose/libraries/RNA_ChunkLibrary.cc` (`setup_base_pair_step_chunks`): a step whose four residues are
  all L-RNA uses `BasePairStepLibrary::mirrored_mini_pose_list()`, the D library poses reflected through
  `MiniPose::reflect_x()`. The chunk is inserted with `copy_dofs`, which takes every DOF from coordinates, so the
  reflection negates all torsions and branch offsets and gives mirrored jumps, as for fragments and jumps.
  A step that mixes D and L residues gets no chunk (one-time warning), and `RNA_DeNovoPoseInitializer` does not
  force a cutpoint inside it. With `-rna:denovo:invert_torsions_for_l_rna false` the old behaviour returns.

### Pucker and ideal geometry
* `core/pose/rna/util.cc`: `assign_pucker()`, `apply_pucker()` (non-idealized branch), and
  `apply_ideal_c2endo_sugar_coords()` (mirrored torsions and sugar coordinates).
* `core/pose/rna/util.cc` `position_cutpoint_phosphate_torsions()`: places the 3' residue's phosphate at a chain
  break with a dummy residue of the same chirality and the mirrored A-form alpha (OP1/OP2 were swapped for L).
* `core/pose/rna/RNA_IdealCoord.{hh,cc}`: L residues (`na_lra` etc.) get mirrored reference
  geometry. Pucker detection for `ANY_PUCKER` is chirality-aware.

### Scoring (every RNA term that depends on handedness)
| Term | Change |
|---|---|
| `rna_torsion` | energy already chirality-aware upstream (mixed D–L bonds score 0, with a warning); its atom derivative (`get_f1_f2`) had the wrong sign of f2 for L-RNA, which made the minimizer strain L sugars. Fixed |
| `rna_sugar_close` | pucker classification uses the mirrored δ |
| syn/anti χ restraints (`add_syn_chi_constraints`, `add_anti_chi_constraints`) | L residues get the mirrored spline table (node q ↔ 37−q) |
| `rna_suite`, `suiteness_bonus` | L:L suites are scored on mirrored torsions, and derivatives are mirrored back. Mixed D:L suites are skipped |
| `rna_jr_suite`, `RNA_SuiteName` | same rules as `rna_suite` |
| `fa_stack`, `stack_elec` | use chirality-normalized base frames. They depend only on cos² of the angle to the base normal, so D:L pairs need no special handling |
| low-res `rna_base_pair`, `rna_base_axis`, `rna_base_stagger`, `rna_base_stack`, `rna_base_stack_axis`, `rna_base_backbone`, and the DMS low-res potential | use chirality-normalized base frames (through `RNA_CentroidInfo::base_stubs()`). For a D:L pair the base-base terms (not `rna_base_backbone` or DMS) use the raw frame for both bases (`base_stub_for_pair()`, below) |
| Base-pair classification (`RNA_BasePairClassifier`) | uses `base_stubs()`, so L:L edges and cis/trans are assigned like their D mirror. D:L pairs use `base_stub_for_pair()` |
| Canonical-pair identity (`possibly_canonical`, `possibly_canonical_strict`, `get_base_pair_atoms` in `core/chemical/rna/util.cc`) | L-RNA AAs are mapped to D (`rna_d_equivalent_aa`), so L:L pairs count as Watson–Crick (N_WC, natWC, rms_stem), helical fragments are picked for L stems, and base-pair constraints get their H-bond atoms. A D:L pair is never canonical and gets no atom-pair constraints |
| `fa_atr`, `fa_rep`, `fa_intra_rep`, `lk_nonpolar`, `geom_sol*`, `hbond_*`, `fa_elec*`, `rna_vdw`, `rna_rg`, `ch_bond` | depend only on distances and angles, so they are mirror-invariant and need no change |
| `ref` | already maps `na_lra` etc. to the D reference energies |

**D:L base pairs.** Bases are planar, so a base's raw frame depends only on its own atoms, not on sugar chirality.
A D:L pair with the same base coordinates as a D:D pair is therefore physically the same base–base arrangement.
Normalizing only the L residue would score it as if one base were flipped over (cis read as trans). So
`core::chemical::rna::base_stub_for_pair()` gives both bases of a D:L pair the raw frame, and the D:D statistics
apply. The base–backbone term and the DMS potential involve the partner's chiral backbone; D:D statistics have no
counterpart for that, and they are left unchanged.

The derivatives of the frame-based terms stay consistent: the normalized frame only replaces `z`
with `−z`, which is still a vector fixed to the base.

### Packing and residue types
* `core/pack/rotamer_set/rna_rotamer_building_functions.cc`: χ rotamer centers, the purine syn
  check, and the 3'-phosphate ε samples are mirrored for L residues.
* `core/chemical/PatchOperation.cc`:
  * the `L_RNA` patch now also mirrors proton-χ (2'-OH) samples, and sets `na_analogue` of each L
    type to its D counterpart (so code that identifies the base through `na_analogue()` works),
  * `ChiralFlipAtoms` also mirrors the ICOORs of the residue connections (LOWER/UPPER) for RNA. They place the
    neighbouring linking atom, and through it the phosphate oxygens and the cutpoint virtual atoms (OVL1/OVL2/OVU1).
  * `SET_ICOOR` negates the dihedral when it patches an L-RNA type. All RNA variant patches are
    written for D-RNA; without this the 5'-terminal phosphate (`LowerRNA`) and other patched atoms
    of L residues were built with D handedness.
* `core/chemical/AA.cc`: `1`–`4` parse to the L-RNA AA types (for `make_pose_from_sequence`), and
  AA names `LRA`, `LRC`, `LRG`, `LUR` exist for params files.
* `core/pose/annotated_sequence.cc` (`residue_types_from_sequence`): an L-RNA residue given by name (as the 1-4
  shorthand produces) at a chain end gets the LowerRNA/UpperRNA terminus variant, like a plain one-letter D residue.
  Other named residue types are taken literally, as before.
* `database/.../patches/nucleic/rna/{VirtualBase,VirtualBaseHeavyAtom,BulgeRNA,Virtual_RNA_Residue_Exclude_Phosphate}.txt`:
  the `AA` selectors also list `LRA LRC LRG LUR`. These patches only change atom types (virtualization), so they
  apply to L unchanged. Without this, rna_denovo stopped on loops it virtualizes (e.g. the UUCG tetraloop).
  Patches that add atoms (BlockStack, Deoxy_O2prime, triphosphate, cyclic phosphate, 7-deaza, protonated
  variants) are still D-only.
* `core/sequence/util.{hh,cc}`: `expand_l_rna_shorthand()`, applied in
  `parse_out_non_standard_residues()` (the `-sequence` / FASTA path used by rna_denovo).

### Pre-existing L-RNA scoring bugs fixed
* `ref`: L-RNA scored 0 because the L AA values lie past the end of the weights list; they now use
  the D weights.
* `rna_torsion`: guanosine-specific χ terms (syn-G bonus, old guanosine χ potential) now apply to L-G.
* Low-res `get_bin()` (RNA_LowResolutionPotential) exited on L-RNA residues.

## Limitations (still D-only or unsupported)

* **Mixed D–L covalent steps.** A D residue bonded to an L residue in the same chain gets no
  fragments across the junction, no `rna_torsion`/`rna_suite` score or `rna_torsion` derivative for that step, and no
  suite name. Put D and L strands in separate chains.
* **Heterochiral (D:L) base pairs or stacks.** No data exists for these. Low-res base-pair and stacking
  terms score them with the normalized frames, which is arbitrary. Jumps between a D and an L residue
  use the unmodified D library. Treat any D:L pairing from these terms as uninformed. Use
  constraints if you know the interface.
* **RNA–protein low-res terms** (`rnp_*`, `HackAro`) are not mirrored. They are D-RNA/L-protein
  statistics, and an L-RNA/L-protein interface has no mirror counterpart in them.
* **`RNA_HelixAssembler`** (`rna_helix`, `-dump_stems`) builds D A-form helices only. For L stems,
  supply your own L helix as a chunk (`-s`). For example, mirror a D helix by negating x in the PDB and renaming
  residues to ` 0A`/` 0C`/` 0G`/` 0U`.
* **Stepwise (SWM/SWA) samplers** and other non-FARFAR RNA samplers that hard-code torsion bins
  were not audited.
* **Chemical shift** terms (`rna_chem_shift`) use cos²-type geometry and are mirror-invariant in
  principle, but this has not been tested.

## Custom L-RNA params

The four custom files (`L_RAD_n.params` etc.) are exact mirrors of the stock `rna_phenix` D params,
so they duplicate the built-in `LRAD`/`LRCY`/`LRGU`/`LURA`. Corrected copies are kept, **unregistered**, in
`database/chemical/residue_type_sets/fa_standard/residue_types/nucleic/rna_l_custom/`. Fixes:
`AA LRA|LRC|LRG|LUR` (was the D code), `NA_ANALOGUE RCY` for C (was `RAD`), and
`IO_STRING ␣0X n` (PDB name ` 0A` etc., matching the built-in types). The one-letter codes stay
`a c g u`: much of the RNA code (jump library, `rna_vdw`, chunk library, base-pair classification)
looks residues up by one-letter code, so a digit code would break it. Use the `1`–`4` shorthand in
sequences instead. Do **not** register them in `residue_types.txt` alongside the `L_RNA` patch: two
types would claim the same AA. To use them,
remove `patches/nucleic/rna/L_RNA.txt` from `patches.txt` first.

## FullModelParameters / FASTA notes

Earlier local changes (bracket handling in `FullModelParameters`, `name1()` remapping in
`check_full_model_info_OK`) are not needed here. Upstream already strips brackets and fills
`non_standard_residue_map` when it parses the FASTA (`core::import_pose::get_sequence_information`,
`core::sequence::parse_out_non_standard_residues`). The shorthand expands to that bracketed form,
and the L types' one-letter codes are `a c g u`, so the `name1()` comparison succeeds. Bracketed
names such as `g[LRGU]` or `3[LRGU]` also work.

## Tests

`test/core/fragment/rna/L_RNA_Chirality.cxxtest.hh` (`fragment/rna` in `core.test.settings`):
* `1234` gives L-RNA types, and the FASTA shorthand expands to the annotated form,
* mirrored torsions give an exact mirror image (all internal distances equal),
* puckers, C2'-endo ideal geometry and `RNA_IdealCoord` agree between D and mirrored L,
* fragment insertion negates torsions for L, leaves D untouched, and skips mixed windows,
* every hi-res (`rna_res_level_energy4` + suite terms) and low-res (`rna_lores`) score term is equal
  for a D pose and its L mirror,
* `rna_torsion` atom derivatives of L match those of the D mirror, also across a D–L junction,
* base-pair-step chunks are mirrored for all-L steps and skipped for mixed steps,
* L·L Watson–Crick pairs count as canonical and heterochiral pairs do not,
* cutpoint variant atoms (OVL1/OVL2/OVU1) and the phosphate after a chain break are exact mirrors,
* named L-RNA residues at chain ends get terminus variants,
* χ restraints (syn/anti) give equal energies for D and the L mirror.
* a D:L pair (an L residue with the same base atoms as the D residue it replaces) gets the same low-res
  base–base energies and the same base-pair classification as the D:D pair (fails without `base_stub_for_pair()`).

Verified: all tests above pass in a release build, as do the existing `chemical`, `sequence`, `pose/rna`
and `test_rna_sfxn` unit tests.

### rna_denovo on a mixed D/L chain
GCAA hairpin with a D stem and an L loop (`-sequence gggc3211gccu -secstruct "((((....))))"`, 200 decoys,
same settings as below). All 200 runs finished, and every decoy has the intended handedness per residue. The
O3'–P bonds at the D→L and L→D junctions have a median of 1.61 Å, the same as elsewhere in the chain (one decoy
in 200 has 1.94 Å at a junction). Stem RMSD to the all-D native (residues 1–4, 9–12) has a median of 1.16 Å,
against 1.11 Å for the all-D hairpin. The median score is −10.2, against −19.7 for all-D, mostly through less
stacking and contact (`fa_stack` +3.6, `fa_atr` +2.6). That is expected: the mirrored loop cannot continue the D
stem's stack, and the junction steps get no fragments, no base-pair-step chunks and no `rna_torsion` term.

### rna_denovo on a two-strand D·L duplex
Szabat 12-mer, D strand `agaaagagaaga` and L strand `424442424424`, parallel, with the 12 cis-WC pairs given as
`-obligate_pair_explicit A:i B:i W W C` (heterochiral pairs cannot be written as dot-bracket). 100 decoys, 3000
cycles. The mechanics work: every decoy has the intended handedness per residue and intact backbone bonds, and
the pairs form (median 12/12; all 12 in 64 decoys). Mixed pairs get the D·D jump unchanged, which is right: the
base and C1' are planar, so the D·D cWW transform gives the D·L cWW base geometry, and the flipped L sugar makes
the strands parallel.

The helix itself is not reproduced. Before the D:L base-pair fix (`base_stub_for_pair()`), mean twist was
about 0° (median −2.5°) against −29° for the left-handed duplex that is stable in MD and +31° for the
right-handed model, and RMSD to those models had medians of 8.7 and 10.0 Å. After the fix the distribution is
much the same (median twist 0°, median RMSD 9.5 Å to either model, all 12 pairs in 71 decoys). The best-scoring
decoys now lean left-handed, though: the top four have negative twist, and the second best is at −27° and 3.2 Å
from the left-handed model.

The full-atom score is not the obstacle. Minimized in Rosetta with the same weights, the left-handed MD helix
keeps its twist (−34°) and all 12 pairs and scores −12.2. That beats the right-handed model (−4.9) and every
decoy (best −6.7). The low-res stage, where the fold is set, still prefers other structures: −152 for the
minimized left-handed helix, −161 for the right-handed model, and down to −168 for untwisted decoys. Each strand
takes its fragments, `rna_torsion` and low-res stacking statistics from homochiral RNA, so the sampling pulls the
D strand toward right-handed and the L strand toward left-handed. Mixed steps also get no base-pair-step chunks.
Sampling, not full-atom scoring, is what has to change for heterochiral helices.

### rna_denovo benchmark, D vs mirrored L
Two hairpins, each folded from sequence as D and as the L mirror (200 decoys each, 3000 cycles,
`-minimize_rna`), with RMSD to the native or its mirror:

| target | | best score | median score | median RMSD (Å) | decoys < 2 Å | RMSD of top-scoring decoy |
|---|---|---|---|---|---|---|
| GCAA hairpin (12 nt) | D | −24.1 | −19.7 | 1.83 | 59% | 1.98 |
| | L | −24.4 | −19.1 | 1.77 | 68% | 1.78 |
| UUCG hairpin (14 nt) | D | −23.7 | −15.7 | 2.73 | 7% | 2.86 |
| | L | −26.9 | −15.5 | 2.68 | 11% | 3.51 |

Two-sample KS tests find no D/L difference in score (p = 0.26, 0.53) or RMSD (p = 0.26, 0.17). Individual
decoys differ because a Monte Carlo step whose energy change is ~1e-13 can be accepted in one handedness and
rejected in the other; with such steps treated as zero, D and L follow the same trajectory.
