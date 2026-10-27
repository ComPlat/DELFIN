# Third-party notices

DELFIN is licensed under LGPL-3.0-or-later: the notice is in [LICENSE](LICENSE),
the licence texts are in [COPYING.LESSER](COPYING.LESSER) and [COPYING](COPYING).
This file lists

1. code by other authors that is part of this repository and of every sdist and
   wheel built from it, with the licence texts its terms require;
2. programs that DELFIN's installers download from their upstream at install
   time, which DELFIN does not distribute;
3. Python dependencies whose licence carries obligations beyond DELFIN's own.

It is an engineering record of what is where and under which licence, not legal
advice. The machine-readable form for every file of the repository is
[REUSE.toml](REUSE.toml) with the licence texts in [LICENSES/](LICENSES).

## 1. Distributed with DELFIN

### 1.1 3Dmol.js 2.5.3

- File: `delfin/dashboard/static/3Dmol-min.js` (minified webpack bundle,
  vendored in commit 2d930031 so that the molecule viewers work without network
  access). Its banner names version 2.5.3; the file is not byte-identical to the
  npm `3dmol@2.5.3` build or to the current `https://3Dmol.org/build/3Dmol-min.js`.
- Upstream: <https://github.com/3dmol/3Dmol.js>, copyright 2014 University of
  Pittsburgh and contributors.
- Licence: BSD-3-Clause for 3Dmol.js. The bundle also contains pako (MIT; its
  zlib port under the Zlib licence), UPNG.js (MIT), iobuffer (MIT), netcdfjs
  (MIT) and an FXAA shader adapted from Rendu (MIT) by way of Mol* (MIT).
  SPDX: `BSD-3-Clause AND MIT AND Zlib`.
- The complete texts of all of these ship beside the bundle in
  `delfin/dashboard/static/3Dmol-min.js.LICENSE.txt`, the file the bundle's
  first line refers to. The 3Dmol.js `LICENSE` (identical in the GitHub tag
  2.5.3 and the npm package 3dmol@2.5.3) reads:

```text
3Dmol.js incorporates code from GLmol, Three.js, and jQuery and
is licensed under a BSD-3-Clause license.

* 3Dmol.js
	Copyright (c) 2014, University of Pittsburgh and contributors
	All rights reserved.

	Redistribution and use in source and binary forms, with or without 
	modification, are permitted provided that the following conditions are met:

	1. Redistributions of source code must retain the above copyright 
	notice, this list of conditions and the following disclaimer.

	2. Redistributions in binary form must reproduce the above copyright 
	notice, this list of conditions and the following disclaimer in the 
	documentation and/or other materials provided with the distribution.

	3. Neither the name of the copyright holder nor the names of its 
	contributors may be used to endorse or promote products derived from 
	this software without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
	"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
	LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
	PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
	HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
	SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
	TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
	PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

* GLmol
 GLmol - Molecular Viewer on WebGL/Javascript (0.47)
  (C) Copyright 2011-2012, biochem_fan
      License: dual license of MIT or LGPL3

  Contributors:
    Robert Hanson for parseXYZ, deferred instantiation



* Three.js 
     https://github.com/mrdoob/three.js

     Copyright (c) 2010-2012 three.js Authors. All rights reserved.

     Permission is hereby granted, free of charge, to any person obtaining a copy
     of this software and associated documentation files (the "Software"), to deal
     in the Software without restriction, including without limitation the rights
     to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
     copies of the Software, and to permit persons to whom the Software is
     furnished to do so, subject to the following conditions:

     The above copyright notice and this permission notice shall be included in
     all copies or substantial portions of the Software.

     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
     IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
     AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
     LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
     THE SOFTWARE.

* jQuery
     http://jquery.org/

     Copyright (c) 2011 John Resig

     Permission is hereby granted, free of charge, to any person obtaining
     a copy of this software and associated documentation files (the
     "Software"), to deal in the Software without restriction, including
     without limitation the rights to use, copy, modify, merge, publish,
     distribute, sublicense, and/or sell copies of the Software, and to
     permit persons to whom the Software is furnished to do so, subject to
     the following conditions:

     The above copyright notice and this permission notice shall be
     included in all copies or substantial portions of the Software.

     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
     LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
     WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
```

### 1.2 OpenTUI patches

- Files: `tui/patches/@opentui%2Fcore@0.6.4.patch` (applied by Bun through
  `patchedDependencies` in `tui/package.json`) and
  `tui/patches/opentui-v0.6.4-pr1665-footer-growth-scrollback.patch` (applied by
  `tui/scripts/build-opentui-native.sh`). Both are shipped in the sdist and,
  through `setup.py`, in the wheel.
- They reproduce source lines of OpenTUI 0.6.4 (<https://github.com/anomalyco/opentui>)
  as diff context and removed lines; the added lines are DELFIN's changes.
- Licence of the OpenTUI parts: MIT. SPDX of the patch files:
  `MIT AND LGPL-3.0-or-later`.

```text
MIT License

Copyright (c) 2025 opentui

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

### 1.3 COD crystal structures in the review sets

- File: `review_sets/kalib_2026_10/structures.xyz.gz` (blinded review set; not
  part of the wheel). Besides MANTA and competitor constructions it holds 600
  crystal structures from the Crystallography Open Database
  (<https://www.crystallography.net/cod/>; S. Gražulis et al., *Nucleic Acids
  Res.* 2012, 40, D420–D427), H positions normalised to standard bond lengths.
- Licence of the COD coordinates: CC0-1.0 (public domain dedication). SPDX of
  the file: `LGPL-3.0-or-later AND CC0-1.0`.
- The COD identifiers are withheld while the set is blinded; they are released
  with the ratings.

## 2. Downloaded at install time, not distributed

DELFIN's installers fetch these programs from their upstream into
`delfin/qm_tools`, `~/.delfin/qm_tools` or `~/.delfin/tools`. They are not part
of any DELFIN sdist or wheel (`scripts/check_distribution.py` fails a build that
contains them). Whoever redistributes a prepared installation that includes them
(a container image, an environment archive, a cluster module) conveys them under
their own licences: for the LGPL/GPL programs that means their licence texts and
their Corresponding Source or a written offer for it.

| Program | Upstream | Licence | Fetched by |
|---|---|---|---|
| xtb | <https://github.com/grimme-lab/xtb> | LGPL-3.0-or-later | `install_qm_tools.sh` (conda-forge) |
| CREST | <https://github.com/crest-lab/crest> | LGPL-3.0-or-later | `install_qm_tools.sh` (conda-forge) |
| xtb4stda, stda 1.6.1 and their parameter files | <https://github.com/grimme-lab/xtb4stda> (release v1.0) | LGPL-3.0-or-later | `install_qm_tools.sh` |
| std2 | <https://github.com/grimme-lab/std2> | LGPL-3.0-or-later | `install_qm_tools.sh` (built from source) |
| g-xTB | <https://github.com/grimme-lab/g-xtb> | GPL-3.0 | `install_qm_tools.sh` (release asset, checked against the release's `.sha256`) |
| DFTB+ | <https://github.com/dftbplus/dftbplus> | LGPL-3.0-or-later | `install_qm_tools.sh` (conda-forge) |
| MOPAC | <https://github.com/openmopac/mopac> | Apache-2.0 | `install_qm_tools.sh` (conda-forge) |
| micromamba | <https://github.com/mamba-org/mamba> | BSD-3-Clause | the installers, when no conda is present |
| Bun | <https://github.com/oven-sh/bun> | MIT; statically links JavaScriptCore/WebKit (LGPL-2.0) and further libraries listed in its `LICENSE.md` | `tui/scripts/install-bun.sh`, `delfin-agent tui` (version- and SHA-256-pinned) |
| npm packages of the terminal window (`@opentui/*`, `solid-js`, `web-tree-sitter`) | npm registry, versions in `tui/bun.lock` | MIT | `bun install --frozen-lockfile` into `~/.delfin/tools` |
| epic-MACE | <https://github.com/EPiCs-group/epic-mace> | GPL-3.0 | `python -m delfin.installer --install epic-mace` (own environment, run as an external program) |
| ORCA | <https://www.faccts.de/> | proprietary (ORCA EULA) | never downloaded; each user installs a licensed copy |

A compiled terminal window (`tui/scripts/build-binary.ts`, `bun build --compile`)
embeds the Bun runtime with JavaScriptCore and every bundled npm package. No
DELFIN release distributes such a binary; one that does needs a notice file
generated by that build (LGPL-2.0 relinking terms of JavaScriptCore, and the
MIT/ISC/BSD notices of the packages).

## 3. Python dependencies with obligations beyond DELFIN's licence

These are declared in `pyproject.toml` and installed by pip from their own
distributions; none of them is part of a DELFIN wheel or sdist. Declaring a
dependency is not distribution, and nothing below restricts a user's own
installation. It matters to whoever redistributes an environment that contains
DELFIN together with them. The job "dependency licences" in
`.github/workflows/packaging.yml` checks the licences of a plain
`pip install delfin-complat` on every pull request against the reviewed
exceptions in `scripts/dependency_licences.toml`.

| Package | Installed by | Licence | What follows |
|---|---|---|---|
| openbabel-wheel (Open Babel 3.1.1) | core | GPL-2.0-only (the wheel's LICENSE is GPL version 2; `openbabel/mol.h` grants version 2 without "or later") | Linked in-process (`from openbabel import pybel`). GPL-2.0-only code cannot be combined with GPL-3.0/LGPL-3.0 code in one distributed work, so a bundle that ships DELFIN together with Open Babel (container image, environment archive, cluster module, installer that carries both) has no set of licence terms that can all be met. Such a bundle should not be distributed; let the target environment install Open Babel itself. |
| PyMuPDF | extra `[pdf]` (also in `[dev]`) | AGPL-3.0, or Artifex commercial licence | Not installed by default. The file browser's PDF viewer, OCR page images and chat previews of PDFs need it; the viewer and the OCR check say so when it is missing, a chat card shows the file without a preview. PDF text, forms and merging work through pypdf without it. AGPL-3.0 section 13 adds a source offer for users of a modified version over a network; whether a dashboard served to other users is such a use is a legal question for whoever operates it. |
| molSimplify | extra `[ai-complex]` | GPL-3.0 | Compatible with LGPL-3.0; a distributed combination that includes it is under GPL-3.0 terms. |
| orca-pi (OPI, FACCTs GmbH) | extra `[orca]` | GPL-3.0 (its PyPI licence text; FACCTs notes that OPI also ships with ORCA under different terms) | Imported in-process by `delfin/opi_witness.py` only when installed. Compatible with LGPL-3.0; a distributed combination that includes it is under GPL-3.0 terms. ORCA itself is proprietary and never bundled. |
| formulas | extras `[agent]`, `[office]` | EUPL-1.1-or-later | The EUPL's compatibility clause lists the GPL; an engineering reading, not reviewed by a lawyer. |
| ase | core | LGPL-2.1-or-later | Compatible. |
| odfpy | extras `[agent]`, `[office]` | Apache-2.0, GPL or LGPL, at the user's choice | Used under Apache-2.0. |

Every other core dependency declares a permissive licence (MIT, BSD, Apache-2.0,
ISC, PSF, MPL-2.0); `stk` and four of its dependencies declare none and are MIT
by their licence files (see `scripts/dependency_licences.toml`).
