/* Editorial content: change text here without touching layout or animation.
   Video examples and editing instructions are in README.md. */
window.DELFIN_CONTENT = {
  hero: {
    eyebrow: 'COMPUTATIONAL CHEMISTRY / AGENTIC SCIENCE',
    title: 'Molecular insight.',
    accent: 'Scientific impact.',
    description: 'From molecular structure to reproducible predictions. DELFIN connects quantum chemistry, molecular design and AI agents in one open research platform.'
  },
  workflow: [
    { title: 'Start with an idea.', label: 'Molecular input', tag: '01 / DEFINE',
      text: 'Bring a SMILES string or a molecular structure. Choose the question you want to investigate, from redox behaviour to excited states.',
      detail: 'Organic molecules & transition-metal complexes', code: 'molecule → research question' },
    { title: 'Give it structure.', label: 'Build & explore', tag: '02 / CONSTRUCT',
      text: 'Generate three-dimensional structures and explore conformers. MANTA brings deterministic coordination construction to metal complexes; established tools cover organic structures and conformer searches.',
      detail: 'MANTA · RDKit · Architector · CREST', code: 'SMILES → 3D structures → conformers' },
    { title: 'Let physics do the talking.', label: 'Calculate', tag: '03 / COMPUTE',
      text: 'Connect quantum-chemistry engines and machine-learning potentials through practical workflows. Agents can configure calculations, follow jobs and help interpret the outputs.',
      detail: 'ORCA · xTB · machine-learning potentials', code: 'structures → calculations → properties' },
    { title: 'Keep the evidence.', label: 'Understand', tag: '04 / LEARN',
      text: 'Inspect results, spectra and method choices. Export reports with provenance, then use your chemical judgement to decide which idea deserves the next calculation.',
      detail: 'Dashboard · DOCX reports · JSON export', code: 'results → evidence → next question' }
  ],
  research: [
    { name: 'Catalysis & energy', number: '01', title: 'Follow the electrons.',
      text: 'Explore redox potentials, spin states and complex stability. Compare molecular candidates to inform research into catalysts and energy-related chemistry.',
      tags: ['Redox potentials', 'Spin states', 'Thermodynamic cycles'], symbol: 'e⁻', note: 'REDOX / STABILITY / REACTIVITY' },
    { name: 'Light & materials', number: '02', title: 'See molecules in a different light.',
      text: 'Investigate absorption, emission and excited-state dynamics. Connect screening and quantum-chemical workflows for photoactive molecules and emissive materials.',
      tags: ['TADF screening', 'Excited-state dynamics', 'UV–Vis spectra'], symbol: 'hν', note: 'ABSORB / EXCITE / EMIT' },
    { name: 'Spectra & structure', number: '03', title: 'Make the invisible readable.',
      text: 'Connect molecular structures to computed spectra. Explore conformer ensembles, vibrational fingerprints and Boltzmann-weighted NMR predictions.',
      tags: ['Ensemble NMR', 'IR spectroscopy', 'Conformer sampling'], symbol: 'δ', note: 'STRUCTURE / ENSEMBLE / SIGNAL' }
  ],
  animals: [
    { name: 'DELFIN', role: 'Connect the science.', description: 'The platform bringing molecular workflows, AI agents and evidence together.', image: 'assets/delfin.png', href: 'https://github.com/ComPlat/DELFIN' },
    { name: 'MANTA', role: 'Build the possibilities.', description: 'Deterministic coordination construction for transition-metal complexes.', image: 'assets/manta.png', href: 'https://github.com/ComPlat/DELFIN/tree/main/delfin/manta' },
    { name: 'ChemDarwin', role: 'Explore what could evolve.', description: 'Enumerate molecular variants. Let the chemist choose what comes next.', image: 'assets/chemdarwin.png', href: 'https://github.com/ComPlat/DELFIN#what-delfin-can-do' },
    { name: 'WEDDELL', role: 'Part of the Animal Loop.', description: 'A growing part of the DELFIN ecosystem.', image: 'assets/weddell.png', href: 'https://github.com/ComPlat/DELFIN' }
  ],
  /* Empty on purpose: the video section appears when real videos are added.
     Use local MP4/WebM or a YouTube embed URL, a local poster, and captions.
     See README.md for the exact format. */
  videos: []
};
