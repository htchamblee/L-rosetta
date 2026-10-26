L-Rosetta: mixed D-RNA / L-RNA modeling
=======================================

A fork of [Rosetta](https://github.com/RosettaCommons/rosetta) that lets **D-RNA and L-RNA (mirror-image
RNA) exist in the same pose** and be folded together by FARFAR (`rna_denovo`).

Stock Rosetta can model an all-L-RNA molecule by mirroring the whole structure, which is enough when the
L molecule is alone. It cannot model a complex that contains both hands at once: every knowledge-based
potential, the fragment library, the jump libraries and the ideal geometries are derived from D-RNA, so
an L residue gets the wrong sign on its torsions, the wrong base frame, and no usable fragments. That
makes the interesting case — a cross-chiral complex, such as an L-RNA aptamer bound to a D-RNA target —
out of reach.

This fork makes each L-RNA residue behave as the exact mirror image of its D counterpart wherever
handedness matters, so a mixed pose scores and folds correctly while a pure D-RNA pose is bit-for-bit
unchanged.

What it adds
------------

* **L-RNA residues from the sequence.** Write them as `1 2 3 4` (L-a, L-c, L-g, L-u) in `-sequence` or
  a FASTA file; `a c g u` are still D-RNA. No custom params files: the L types are generated from the D
  types by the `L_RNA` patch, with PDB names ` 0A`, ` 0C`, ` 0G`, ` 0U`.
* **Chirality-aware fragment insertion, jumps and pucker.** Vall fragments are sign-inverted into L
  windows; base-pair jumps are conjugated through the reflection; ideal geometries are mirrored.
* **Chirality-aware scoring.** Every RNA term that depends on handedness — base pairing and stacking,
  torsions and their derivatives, the low-resolution potentials — reduces to the original arithmetic on
  a pure D pose.
* **Heterochiral base-pair geometry.** A D·L Watson-Crick pair is scored and classified in the raw base
  frame rather than the chirality-normalized one, which is what makes a cis pair read as cis. An
  all-*anti* D·L cis-WC pair requires **parallel** strands, unlike a conventional D·D pair.
* **A parallel heterochiral base-pair-step library**, in
  `database/sampling/rna/base_pair_steps/canonical/heterochiral_parallel/`, built from MD of a
  left-handed parallel D·L duplex. Stock Rosetta detects antiparallel steps only, so a parallel D/L
  duplex previously got no step chunks at all.
* **Terminus and virtualization variants** for named L-RNA residues, so chain ends, bulges and virtual
  bases work on both hands.

Roughly 2,400 lines across 57 files, on top of the Rosetta mainline. **Full detail, including what is
still D-only, is in
[`source/src/core/fragment/rna/L_RNA_SUPPORT.md`](source/src/core/fragment/rna/L_RNA_SUPPORT.md)** —
read that before changing anything here.

Building
--------

Rosetta is built with the bundled scons; nothing needs to be installed system-wide. You need a C++17
compiler, Python 3, and zlib. No build extras are required — this is a plain `default` build.

```sh
cd source
./scons.py -j8 mode=release bin/rna_denovo.default.linuxgccrelease
```

* scons fetches the external submodules it needs; you do not have to `git submodule update` by hand.
* `-j8` on a machine with 32 GB is comfortable. Each compile job can want over a GB, so do not set `-j`
  much above the core count.
* Expect an hour or more for the first build, and about 1.3 GB of build output. The repository including
  the 2.3 GB `database/` comes to roughly 8 GB checked out and built.
* Drop the `bin/...` target to build every application, which takes considerably longer and is not
  needed for RNA work.
* The binary lands in `source/bin/` as a symlink; the `database/` directory must stay alongside it, or
  pass `-database <path>` explicitly.

Running
-------

A two-strand D·L duplex, D strand first:

```sh
source/bin/rna_denovo.default.linuxgccrelease \
    -fasta duplex.fasta -secstruct_file duplex.secstruct \
    -nstruct 100 -cycles 3000 -minimize_rna true
```

with a FASTA whose L strand uses digits (`1 2 3 4` = L-a, L-c, L-g, L-u, so the L strand below reads
`cuggga` paired against the D strand's `aacccg`):

```
>duplex A:1-6
aacccg
>duplex B:7-12
243331
```

Everything chirality-related turns on by itself once the pose contains L-RNA. The one override,
`-rna:denovo:invert_torsions_for_l_rna false`, restores the old behaviour and is only correct for a pure
L-RNA pose with an already-mirrored vall.

Tests
-----

```sh
cd source
./scons.py -j8 mode=release cat=test
test/run.py --database ../database -1 core:L_RNA_ChiralityTests
```

`source/test/core/fragment/rna/L_RNA_Chirality.cxxtest.hh` covers the mirroring invariants: a D pose and
its reflection must score identically, torsions must negate, and D·L pairs must classify in the right
geometry.

Upstream Rosetta
----------------

Rosetta is maintained by the RosettaCommons and is **not** open source — commercial use requires a
separate licence from the University of Washington. See [LICENSE.md](LICENSE.md), and
<https://www.rosettacommons.org> for the project itself. The upstream repository is
<https://github.com/RosettaCommons/rosetta>; its documentation is at
<https://docs.rosettacommons.org>, and build options are described under
`source/tools/build/` (`site.settings`, `user.settings`).

If you publish work that uses Rosetta, cite it as described in
[CITING_ROSETTA.md](CITING_ROSETTA.md).
