# Parallel heterochiral (D/L) base-pair steps

Steps of a parallel D-RNA/L-RNA duplex with all-anti cis Watson-Crick pairs. Each entry has four residues, in
the order (i, i+1, j, j+1): the D strand first (i, i+1), then the L strand (j, j+1). The pairs are i.j and
(i+1).(j+1). The file name is the D-strand dinucleotide, then the L-strand dinucleotide. A step whose first
strand is L gets the mirror image of an entry (`BasePairStepLibrary::mirrored_mini_pose_list`).

Source: Amber MD (OL3 + CP corrections, Dudek L-RNA parameters), 300 ns, of the Szabat 12-mer duplex
D-`AGAAAGAGAAGA` / L-`UCUUUCUCUUCU`, built as a left-handed cis-WC helix (par-DL-cWW-left, run 1). Frames at
1 ns. The interior steps (pairs 2-11) gave GA/CU, AA/UU and AG/UC. Steps were kept when both pairs had a WC
N1/N3 distance under 3.4 A and a C1'-C1' distance of 9.5-11.5 A, and every second kept frame was used.
The mirror image of each frame (strands swapped so the D strand stays first) gave CU/GA, UU/AA and UC/AG.
Those steps are right-handed, as the mirror of a left-handed helix must be.

Fold tree of each entry: EDGE 1 2, jump 1 -> 3 (the first pair, default RNA jump atoms), EDGE 3 4.
Sequences not listed here have no heterochiral steps; rna_denovo then samples those steps with fragments only.

G·C entries (`cc_gg`, `ac_ug`)
------------------------------

Source: Amber MD (OL3 + CP, Dudek L-RNA parameters, microMg Mg2+, TIP3P), two independent 300 ns
replicas of a parallel D·L interface duplex — D-`AACCCG` against L-`CUGGGA`, pairs 2·8 (A·U) and
3·9 / 4·10 / 5·11 (C·G), with residues 1, 6, 7 and 12 unpaired. Frames at 500 ps, both replicas pooled,
no stride. Same cis-WC filter as the entries above (WC N1/N3 under 3.4 A, C1'-C1' 9.5-11.5 A).

| entry | conformers | note |
|---|---|---|
| `cc_gg` | 2387 | from pairs 3·9/4·10 and 4·10/5·11, formed 99.2 % and 99.7 % of frames in the two replicas |
| `ac_ug` | 212 | needs the A·U pair 2·8, formed 33.7 % in one replica and 1.7 % in the other |

These are the first G·C steps in this library; the entries above are all A/U.

Two caveats on the counts. `cc_gg` is about five times the density of the A/U entries, and its frames
come from a pair that was essentially always formed, so they sample one basin densely rather than
covering more geometry. `ac_ug` is thinner than it looks: all 212 come from the single replica in which
the A·U pair was formed at all, and within that run they fall in a couple of long formed blocks, so
their diversity is lower than the count suggests. The parallel heterochiral A·U pair is weak and 300 ns
per replica does not determine how weak.

`gg_cc` (2387) and `ug_ac` (212) are the reciprocals, made from the two entries above by reflecting each
step and swapping its strands so the D strand is first again (`steplib/mirror_steps.py`). They are
right-handed, as the mirror of a left-handed step must be, and a construct whose L strand comes first in
the pose asks for them rather than for `cc_gg` / `ac_ug`.

Note for anyone regenerating these: BinarySilentStruct stores a one-letter sequence, and L residues
share the letters `a c g u` with D, so chirality does NOT survive the roundtrip out of a silent file -
every residue comes back labelled as D. The coordinates are intact, so the mirror script measures
chirality from the signed volume at C4' rather than trusting the residue names.
