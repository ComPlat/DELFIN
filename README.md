# DELFIN

> **Prereqs**
>
> * ORCA **6.1.0** in your `PATH` (`orca` and `orca_pltvib`)
> * Optional: `crest` (for CREST workflow), **xTB** if used (`xtb` and `crest` in `PATH`)
> * Python **3.9+** recommended

---

## Install

From the `delfin` folder (the one containing `pyproject.toml`):


editable (dev) install
```bash
pip install -e .
```
or regular install
```bash
pip install .
```
global install 
```bash
sudo python3 -m pip install . --break-system-packages
```

This exposes the console command **`delfin`** and enables `python -m delfin`.

---

## Quick start

Create a working directory with at least these two files:

* `CONTROL.txt` — your control/config file
* `input.txt` — the starting geometry (XYZ body without the first two header lines)
* starting from a `XYZ` file is optional

Then run:

from the directory that contains CONTROL.txt and input.txt
```bash
delfin
```
alternatively
```bash
python -m delfin
```

**Convenience (`--define`)**

create CONTROL.txt and an empty input.txt, then exit
```bash
delfin --define
```
convert an existing XYZ → input.txt (drops the first two header lines) and write CONTROL.txt, then exit
```bash
delfin --define=your.xyz
```
overwrite existing CONTROL.txt / input file when defining
```bash
delfin --define=mycoords.txt --overwrite
```
clean up intermediates from previous runs and exit
```bash
delfin --cleanup
```
show options and prerequisites
```bash
delfin --help
```
Re-parse existing outputs and (re)run only external jobs with missing/incomplete .out files
```bash
delfin --recalc
```
The tool writes results and reports into the current working directory,
e.g. `DELFIN.txt`, `OCCUPIER.txt`, and step folders.
---

## Project layout

```
delfin/
  __init__.py
  __main__.py       # enables `python -m delfin`
  cli.py            # entry point orchestrating the full workflow & CLI flags (--define, --cleanup, ...)
  main.py           # optional small loader (may delegate to cli.main)
  define.py         # CONTROL template generator (+ .xyz → input.txt conversion)
  cleanup.py        # delete temporary files
  config.py         # CONTROL.txt parsing & helpers
  utils.py          # common helpers (transition metal scan, basis-set selection, electron counts)
  orca.py           # ORCA executable discovery & runs
  imag.py           # IMAG workflow (plotvib helpers, imaginary-mode loop)
  xyz_io.py         # XYZ/ORCA-input read/write helpers
  xtb_crest.py      # xTB / GOAT / CREST / ALPB solvation workflows
  energies.py       # extractors for energies (FSPE, Gibbs, ZPE, electronic energies)
  report.py         # writing DELFIN.txt and OCCUPIER.txt
  occupier.py       # OCCUPIER workflow (sequence execution + summary)
  copy_helpers.py   # file passing between OCCUPIER steps (prepare/copy/select)
```
---
## Typical workflow switches (in CONTROL.txt)

* `method = OCCUPIER | classic | manually`
* `calc_initial = yes | no`
* `oxidation_steps = 1,2,3` (string; steps to compute)
* `reduction_steps = 1,2,3` (string; steps to compute)
* `E_00 = yes | no`
* `absorption_spec = yes | no`
* xTB/CREST flags: `XTB_OPT`, `XTB_GOAT`, `CREST`, `XTB_SOLVATOR`
* Basis/functional/solvent/RI flags (as before)

---

## Troubleshooting

* **`CONTROL.txt` not found**
  DELFIN exits gracefully and tells you what to do. Create it via `delfin --define` (or copy your own).

* **Input file not found**
  DELFIN exits gracefully and explains how to create/convert it.
  If you have a full `.xyz`, run: `delfin --define=your.xyz` → creates `input.txt` (drops the first two header lines) and sets `input_file=input.txt` in CONTROL.

* **ORCA not found**
  Ensure `orca` is callable in your shell: `which orca` (Linux/macOS) or `where orca` (Windows).
  Add the ORCA bin directory to your `PATH`.

* **`ModuleNotFoundError` for internal modules**
  Reinstall the package after copying files:


* **CREST/xTB tools missing**
  Disable the corresponding flags in `CONTROL.txt` or install the tools and put them in `PATH`.

---

## Dev notes

* Update CLI entry point via `pyproject.toml`
  `"[project.scripts] delfin = \"delfin.cli:main\""`
* Build a wheel: `pip wheel .` (inside `delfin/`).
* Run tests/workflow locally using a fresh virtual environment to catch missing deps.

---
## References

The generic references for ORCA, xTB and CREST are:

- Frank Neese. The ORCA program system. *Wiley Interdiscip. Rev. Comput. Mol. Sci.*, 2(1):73–78, 2012. doi:<https://doi.wiley.com/10.1002/wcms.81>.
- Frank Neese. Software update: the ORCA program system, version 4.0. *Wiley Interdiscip. Rev. Comput. Mol. Sci.*, 8(1):e1327, 2018. doi:<https://doi.wiley.com/10.1002/wcms.1327>.
- Frank Neese, Frank Wennmohs, Ute Becker, and Christoph Riplinger. The ORCA quantum chemistry program package. *J. Chem. Phys.*, 152(22):224108, 2020. doi:<https://aip.scitation.org/doi/10.1063/5.0004608>.
- Christoph Bannwarth, Erik Caldeweyher, Sebastian Ehlert, Andreas Hansen, Philipp Pracht, Jan Seibert, Sebastian Spicher, and Stefan Grimme. Extended tight-binding quantum chemistry methods. *WIREs Comput. Mol. Sci.*, 11:e1493, 2021. doi:<https://doi.org/10.1002/wcms.1493>. *(xTB & GFN methods)*
- Philipp Pracht, Stefan Grimme, Christoph Bannwarth, Florian Bohle, Sebastian Ehlert, Gunnar Feldmann, Jan Gorges, Max Müller, Timo Neudecker, Christoph Plett, Sebastian Spicher, Pascal Steinbach, Piotr A. Wesołowski, and Fabian Zeller. CREST — A program for the exploration of low-energy molecular chemical space. *J. Chem. Phys.*, 160:114110, 2024. doi:<https://doi.org/10.1063/5.0197592>. *(CREST)*

Please always check the output files—at the end, you will find a list of relevant papers for the calculations. Kindly cite them. Please do not only cite the above generic references, but also cite in addition the
[original papers](https://www.faccts.de/docs/orca/6.0/manual/contents/public.html) that report the development and ORCA implementation of the methods DELFIN has used! The publications that describe the functionality implemented in ORCA are
given in the manual.



# Dependencies
 
To use DELFIN, you must be authorized to use ORCA 6.1.0. You can download the latest version of ORCA here:  
https://orcaforum.kofo.mpg.de/app.php/portal  

***ORCA 6.0.1 is free for academic use. Please note the license conditions of ORCA when using it.***  
https://www.faccts.de/

***xTB is free for academic use under the GNU General Public License (GPLv3).***  
The code and license information are available here: https://github.com/grimme-lab/xtb

***CREST is free for academic use under the GNU General Public License (GPLv3).***  
The code and license information are available here: https://github.com/crest-lab/crest


---

## Please cite

If you use DELFIN in a scientific publication, please cite:

Hartmann, M. et al. (2025). *<TITLE OF THE PAPER>*. <JOURNAL>, <VOLUME>(<ISSUE>), <FIRST–LAST>. https://doi.org/<DOI>

### BibTeX
```bibtex
@article{<your-key>,
  author  = {Hartmann, Maximilian and <Coauthors>},
  title   = {<TITLE OF THE PAPER>},
  journal = {<JOURNAL>},
  year    = {2025},
  volume  = {<VOLUME>},
  number  = {<ISSUE>},
  pages   = {<FIRST--LAST>},
  doi     = {<DOI>},
  url     = {https://doi.org/<DOI>}
}
```
## License

This project is licensed under the GNU Lesser General Public License v3.0 or later (LGPL-3.0-or-later).

You should have received a copy of the GNU Lesser General Public License along with this repository in the files `COPYING` and `COPYING.LESSER`.  
If not, see <https://www.gnu.org/licenses/>.

Non-binding citation request:  
If you use this software in research, please cite the associated paper (see [CITATION.cff](./CITATION.cff)).


  