QM Tools Bundle
Date: 2026-03-11

Purpose
- Provide one local entry point for the QM stack needed by DELFIN and ChemDarwin.
- Keep everything accessible from a single PATH prefix under:
  delfin/qm_tools/bin
- Keep installer scripts and notes in git, but leave local binaries/build products untracked.

What is inside
- Native local binaries downloaded into qm_tools/bin:
  - xtb4stda
  - stda_v1.6.1
  - stda -> stda_v1.6.1
- Local entry points to already verified system tools:
  - xtb
  - crest
  - std2
  - dftb+

Why some tools are symlinked instead of copied
- CREST, xTB, std2, and DFTB+ are already installed and working on this server.
- Repackaging those binaries here adds no value and makes updates harder to track.
- The local bin directory still gives DELFIN one stable tool location to target.

Activation
- source delfin/qm_tools/env.sh

Installer
- delfin/qm_tools/install_qm_tools.sh
- The installer downloads the freely available tools and installs xtb/crest/dftb+ into a local micromamba/conda environment if needed.
- ORCA is intentionally not managed in qm_tools.
- Examples:
  - reuse system tools when available:
    USE_SYSTEM_TOOLS=1 delfin/qm_tools/install_qm_tools.sh
  - force a local std2 source build:
    INSTALL_STD2_FROM_SOURCE=1 delfin/qm_tools/install_qm_tools.sh

Health check
- delfin/qm_tools/check_qm_tools.sh

Relevant notes
- DELFIN integration notes:
  delfin/qm_tools/DELFIN_INTEGRATION_NOTES.txt
