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
