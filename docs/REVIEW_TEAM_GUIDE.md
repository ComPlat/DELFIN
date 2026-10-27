# Rating structures as a team: a guide for reviewers

Thank you for helping. You will look at computer-built 3D structures of
molecules (mostly metal complexes) one at a time and say for each whether a
chemist could accept it (**Pass**) or whether something is wrong
(**Block**), and if so what and where. Your judgement is used to measure how
well an automatic checker sees defects, and to find what the structure
builder should fix next.

Anyone can help, without asking first: **just pick your name**. Several
people share the work: each of you gets a block of 500 structures of your own,
plus 50 that *everyone* rates (so we can measure how much reviewers agree).
You do not need to know which ones those are. Each rating is saved the moment
you give it, and you can stop and continue at any time.

A video walk-through is in [review_team_guide.webm](review_team_guide.webm).

![Rating one structure](review_team_guide.gif)

## 1. Get the set (once)

You need a clone of DELFIN and the dashboard installed (see the
[README](../README.md#-quick-start)). Then unpack the structures of the set
into a folder of your own:

```bash
cd DELFIN                                    # your clone
git pull                                     # the newest assignment
delfin review unpack review_sets/kalib_2026_10 ~/review_kalib_2026_10 --name <your name>
```

`~/review_kalib_2026_10` now holds one `.xyz` file per structure, the
category list (`kategorien.txt`) and the shares given out so far
(`assignment.json`). If your name is listed there, that is your share.
Otherwise a block of 500 structures nobody has yet is drawn for you, from
your name (the same name always gets the same block), and kept in your own
`review_<your name>.json`. Pick a name nobody else in the project uses, and
keep it. Running `unpack` again never touches your ratings. (Without
`--name` the block is drawn the first time you open the folder in the
Review tab.)

## 2. Open your share in the Review tab

1. Start the dashboard: `delfin-voila` (it prints a link with an access
   token; open it in the browser) and go to the tab **Review**.
2. Under **Open a review**: put the folder (`~/review_kalib_2026_10`) into
   **Open**, and your **Name** (always the same name when you come back).
3. **Categories** fills itself from `kategorien.txt` when the folder is
   opened: 18 entries, `H` first and `SYMMETRIE` last. Leave **Blinded**
   unticked.
4. Press **Open**. The status line says *Your share: N structures*; the
   counter on the right, *x / N rated*, the file list and **Next unrated**
   cover only your share.

## 3. Rate a structure

Turn and zoom the molecule in the viewer (mouse), then:

* **P** (or the green button) = **Pass**: the molecule could look like this
  in reality; you would accept it as a starting geometry for a calculation
  or as a figure.
* **B** (or the red button) = **Block**: at least one defect a chemist would
  not accept.

For **Block**, first choose the categories (several are allowed, see the
table below), and mark where each defect is:

1. Press the key of a category (or click it). It is now *active* (dark
   frame).
2. Click the one or two atoms where this defect sits. They get a ring in the
   colour of the category; click again to remove.
3. Choose the next category and click its atoms, and so on. The key of a
   category already chosen but not active makes it active again; the key of
   the active category removes it (with its atoms). The x button removes all
   marks.
4. Press **B**.

**Mark every visible defect, not only the worst one.** Otherwise we cannot
measure, category by category, what the checker finds and what it misses.

A different isomer or conformer than you would expect is **not** a defect,
as long as the structure itself is chemically sound.

### Keys

| Key | Does |
|---|---|
| **P** / **B** | Pass / Block (saves and moves to the next unrated structure) |
| **1**-**9**, **0**, **Shift+1** ... **Shift+8** | toggle the categories, in the order of the table |
| **U** | unsure: writes `unsicher` into the comment |
| **N** | next unrated structure |
| **Left** / **Right** | previous / next structure |
| **H** | hide / show the hydrogen atoms |
| **R** | reset the view |
| **Esc** | leave the comment field (the keys work again) |

The field above the file list jumps to a structure (`T0123` or `123`, then
Enter). The list marks each structure with a tick (pass), a cross (block) or
a dot (open).

### Not sure?

Press **U** (it writes `unsicher` into the comment) and still decide **P** or
**B** as well as you can. A free comment is always welcome, and for the
category `SONST` (other) it is required: say what you saw, new categories
come from those comments.

## 4. The 18 categories

| Key | Code | Meaning |
|---|---|---|
| 1 | `H` | a hydrogen is missing, stretched, at a wrong angle, points into the metal or sits on the metal |
| 2 | `LAENGE_ML` | a metal-ligand or metal-metal distance visibly too long or too short (the ligand is still bound) |
| 3 | `DURCHDRINGUNG` | interpenetration: ligands or rings stuck through each other, a bond threaded through a ring |
| 4 | `CLASH` | atoms too close without interpenetration (also H...H) |
| 5 | `METALLKONTAKT` | a heavy atom that is not a donor sits at the metal (a hydrogen at the metal: `H`) |
| 6 | `TOPO` | topology: a bond missing or one too many, a ligand clearly torn off, a wrong donor atom, a fragment or counter-ion floating away or stuck inside another |
| 7 | `HYB` | hybridisation / angle: a wrong angle, sp3 flat, sp2 pyramidal, a linear group bent (also one angle at a donor) |
| 8 | `PLANAR` | an aromatic ring or sp2 framework not flat, an amide, ester or conjugated bond twisted |
| 9 | `POLYEDER` | the coordination geometry at the metal is implausible |
| 0 | `SONST` | something else: **always write a comment** |
| Shift+1 | `DONOR` | a whole donor turned away from the metal, its lone pair not pointing at the metal, the metal out of the donor plane |
| Shift+2 | `LAENGE_KOV` | a covalent bond too long or too short |
| Shift+3 | `KOLLAPS` | collapse: a ring or ligand balled up, the molecule too compact |
| Shift+4 | `HAPTO` | a pi ligand slipped, wrong hapticity |
| Shift+5 | `STEREO` | wrong configuration (R/S, E/Z, atropisomer) |
| Shift+6 | `KONFORMER` | an unrealistic conformation (eclipsed, a ring bent, a boat without reason) |
| Shift+7 | `FORMEL` | formula: atoms missing or too many, an H on the wrong atom (tautomer, protonation) |
| Shift+8 | `SYMMETRIE` | chemically equal parts visibly built unequal (equal ligands at different distances, a Cp ring tilted) |

Where two categories are close:

* stuck through each other = `DURCHDRINGUNG`; only too close = `CLASH`;
* one angle at a donor = `HYB`; the whole donor turned away = `DONOR`;
* still bound but too long = `LAENGE_ML`; clearly gone = `TOPO`;
* a hydrogen at the metal = `H`; a heavy atom at the metal = `METALLKONTAKT`.

## 5. Saving, stopping, continuing

Every rating is saved at once into `review_<your name>.json` in your folder
(the line *saved ...* under the comment shows it). Stop whenever you like:
open the same folder with the same name again and the tab continues at the
first structure you have not rated yet. You do not have to finish your
share; every finished part is useful.

## 6. Send your ratings back

From your clone, in a terminal:

```bash
delfin review submit ~/review_kalib_2026_10 --name <your name>
```

This copies your `review_<your name>.json` to
`ratings/kalib_2026_10/review_<your name>.json` in a commit of its own on a
new branch and opens a pull request with the GitHub command line `gh` (into
your fork of DELFIN when you cannot push to the project; `gh auth login`
once beforehand). The pull request adds only your own file, so reviewers
never overwrite each other. You can submit again later; the new file
replaces your earlier one.

Without `gh` or GitHub access: the command says what is left to do, or
simply send the file `review_<your name>.json` to the coordinator.

## For the coordinator

```bash
delfin review assign review_sets/kalib_2026_10 Anna Ben --ratings ratings/kalib_2026_10   # reserve blocks (optional), then commit
delfin review merge  review_sets/kalib_2026_10 --ratings ratings/kalib_2026_10            # CSV + agreement
delfin review merge  review_sets/kalib_2026_10 --ratings ratings/kalib_2026_10 --split calibration
```

Nobody has to be added before helping: a name not in `assignment.json`
draws its own block (`self_block`, 500) in its review file, so helpers'
pull requests never touch a common file. `assign` is only needed to reserve
blocks for named reviewers: it keeps every share already given out and gives
new names only structures nobody has yet (`--block`, default 500), and it
leaves out the self-drawn blocks it finds in `--ratings`. The shared
structures (`--overlap`, default 50) and the test part are fixed when the
assignment is first drawn. Two self-drawn blocks may overlap; those frames
are simply rated twice.

`merge` writes one CSV with every rating of every reviewer (columns `rater`,
`shared`, `share` = shared / own / extra, `split`) and prints Cohen's kappa
(two reviewers) or Fleiss' kappa (more; the number of reviewers may differ
per frame) with the share of frames on which all agree and the number of
frames, for ok/defect and for every category, on every frame rated by two or
more. Ratings of the **test** part are a pre-registered hold-out: never use
them to tune anything.
