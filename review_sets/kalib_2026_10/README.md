# Review set kalib_2026_10

4996 three-dimensional structures of molecules, mostly transition-metal
complexes, to be rated by eye in the Review tab of the DELFIN dashboard. The
ratings calibrate the automatic structure checker and show which defects the
structure builder should fix next. How to take part:
[docs/REVIEW_TEAM_GUIDE.md](../../docs/REVIEW_TEAM_GUIDE.md).

* `structures.xyz.gz` -- all structures, one XYZ frame each; the comment line
  is only the structure's id (`T0001` ...). `delfin review unpack` writes
  them as single files.
* `kategorien.txt` -- the 18 defect categories the Review tab offers.
* `assignment.json` -- the 50 structures everyone rates, the test part, and
  blocks reserved for named reviewers (`delfin review assign`).

**Just pick your name.** Anyone can help without being added:
`delfin review unpack review_sets/kalib_2026_10 <folder> --name <your name>`
(or the first opening in the Review tab) draws a block of 500 structures
nobody is assigned, from your name, plus the 50 everyone rates. The block is
kept in your own `review_<name>.json`, so contributions never edit a common
file.

**What is in it.** Structures built by MANTA (several program versions)
mixed with experimental crystal structures from the Crystallography Open
Database (COD, https://www.crystallography.net, released under CC0), in a
random order. Some structures appear twice, rotated differently. The
selection is deliberately broad (metals, coordination numbers, ligand types,
and every kind of construction defect as well as none). Every structure is
centred, randomly rotated and its atoms are sorted the same way; hydrogens
of the crystal structures are set to the usual neutron bond lengths, so the
source cannot be read off the file. Where each structure comes from is kept
back until the ratings are in.

**Calibration and test.** Before any rating, a fixed fifth of the structures
(`split.test` in `assignment.json`, seed `20261027:split`; a structure shown
twice stays in one part) was set aside as the test part. Ratings of the test
part are never used to tune anything; they only measure the tuned result.

Ratings come back as `ratings/kalib_2026_10/review_<name>.json`
(`delfin review submit`).
