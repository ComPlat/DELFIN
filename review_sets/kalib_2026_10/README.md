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
* `assignment.json` -- the 50 structures everyone rates, the test part,
  blocks reserved for named reviewers (`delfin review assign`), the
  parameters of the draw (seed, block of 500, redundancy 0.15) and the set
  version.

**Just pick your name.** Anyone can help without being added:
`git pull`, then
`delfin review unpack review_sets/kalib_2026_10 <folder> --name <your name>`
(or the first opening in the Review tab) draws a block of 500 structures
plus the 50 everyone rates. Structures nobody has rated or reserved come
first; about 15 % of every block are structures meant for two or three
reviewers, spread over the set, so that every reviewer's ratings overlap
with many others'. The block is kept in your own `review_<name>.json`, so
contributions never edit a common file. `delfin review claim <folder> --name
<your name>` reserves it for 28 days (a pull request that adds only
`ratings/kalib_2026_10/claims/<name>.json`). `delfin review status
review_sets/kalib_2026_10` shows how far the set is.

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

*Change of 27.10.2026, before any rating was sent back:* the original split
kept every structure shown twice in one part, but not the structures that
are geometrically identical for another reason (the same molecule built
alike under different entries). 14 groups of such identical structures
straddled the two parts; each was moved wholly into the part that held most
of it (ties decided by the test size and the seed). 14 structures moved: 6
into the test part, 8 into calibration; the test part now holds 997 of the
4996 structures. Nothing else changed: the 50 structures everyone rates and
the named block are as before (`split.changes` in `assignment.json`).

Ratings come back as `ratings/kalib_2026_10/review_<name>.json`
(`delfin review submit`, optionally with `--credit "<name to be listed>"`),
reservations as `ratings/kalib_2026_10/claims/<name>.json`. Every such pull
request is checked automatically (`.github/workflows/ratings.yml`).
