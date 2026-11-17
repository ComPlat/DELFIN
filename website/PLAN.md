# DELFIN website — review draft

## Direction

An English research-platform website with an ambitious environmental mission.
Headline proposal: “Small molecules. Planet-sized purpose.” All editorial copy
is a proposal for discussion, not a new official statement of project policy.

Visual identity comes from the actual DELFIN dolphin logo, dashboard wordmark
(`#628FBF`), terminal theme (`#6290C0`, `#8FB2DA`) and pale blue login surfaces.
Use generous typography, alternating light and deep-blue chapters, thin scientific
linework, a procedural molecular illustration, ocean-like waves, scroll reveals,
and restrained pointer interaction. Respect reduced-motion preferences.

## Narrative

1. Purpose and an immediate explanation of the platform.
2. Interactive journey: molecular input → structure → computation → evidence.
3. Research applications: catalysis/energy, photophysics, spectroscopy.
4. Agent, dashboard and CLI: different ways into the same research platform.
5. Mission: science for humanity and the planet, grounded in the agent principles.
6. Optional growing video library, configured separately from the layout.
7. Installation, documentation, contribution and citation.

## Evidence behind the draft

- `README.md`: platform scope, workflows, interfaces, installation, limitations.
- `delfin/agent/pack/shared/principles_addendum.md`: mission and operating principles.
- `tui/src/theme/theme.ts`: restrained steel-blue palette.
- `delfin/dashboard/__init__.py`: original wordmark colour.
- `delfin/voila_assets/login_templates/_delfin_card_style.html`: light/dark surfaces.
- `delfin/logo/DELFIN_logo.png`: unmodified original brand asset.
- `https://github.com/hmaximilian/hmaximilian/blob/main/README.md`: Max's
  personal voice, “Sense → Think → Act → Learn”, and “Stay in the Animal Loop”.
  The profile presents DELFIN, MANTA, WEDDELL and ChemDarwin together. The draft
  uses their original local logos and gives this family its own animated chapter.

## Feedback incorporated during the first draft

- Strong personal identity, not an anonymous software landing page.
- Rich motion: rotating molecular illustration, orbital paths, ocean-like waves,
  moving loop, floating animal marks, scroll reveals and interactive chapters.
- Global motion pause plus reduced-motion system preference; no autoplay video.
- Content remains provisional so Max can add, remove and reword it.
- A growing media library: configurable videos, local posters, captions and
  click-to-load playback. Empty until real recordings are supplied.

No invented performance statistics, testimonials, scientific results or claims
of environmental impact. Molecular art and interface previews are explicitly
illustrative. ChemDarwin's selection loop remains chemist-led.

## Second review: scientific presentation

The initial illustrative molecular artwork was replaced at Max's request with
actual MANTA-generated fac-Ir(ppy)3 XYZ coordinates. The structure explorer has
61 atoms, element colours, pointer/keyboard rotation, a hydrogen toggle, reset,
and a downloadable XYZ. `assets/STRUCTURE.md` records provenance and limitations.
The source geometry is not DFT-optimised; no property predictions are implied.

The design now uses a dark technical hero, clearer sans-serif typography, quieter
motion and a static scientific loop band instead of a scrolling marketing strip.
The proposed hero is “Molecular insight. Scientific impact.” The personal “m.”
signature and all personal profile links/teasers were removed. The logo reads
“DELFIN” without a trailing full stop. Earlier concept notes above record the
initial proposal, not the final visual implementation.

## Growth and publication

Buildless static HTML/CSS/JavaScript, compatible with a GitHub Pages project path.
Core content lives in `content.js`; videos are an extensible array there. New
chapters can be added as HTML sections. Assets are local; no tracking, external
fonts, remote scripts or third-party video requests on initial page load.

Review locally first. No commit, push, Pages configuration or deployment is part
of this draft. Once approved, publish only the website directory using GitHub
Pages; keep internal repository documents out of the published artifact.
