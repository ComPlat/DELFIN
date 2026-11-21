"""Portable offline 3Dmol view plus explicit PyMOL/ChimeraX entry scripts."""
from __future__ import annotations

import json
from pathlib import Path
import shutil

from delfin.tools.bio import BindingBox, read_structure
from delfin.tools.bio.artifacts import file_hash, read_json, write_json
from delfin.tools.bio.processing import load_ligands
from delfin.tools.bio.structures import select_structure


_PAGE = '''<!doctype html>
<html lang="en"><meta charset="utf-8"><title>DELFIN molecular view</title>
<style>body{font:16px system-ui;margin:0;background:#f5f4ee;color:#25352f}header{padding:18px;display:flex;gap:20px;align-items:center;flex-wrap:wrap}#viewer{height:78vh;position:relative}#error{color:#a22;padding:1em}label{display:flex;gap:8px}footer{padding:12px}</style>
<header><strong>DELFIN molecular view</strong><label>Chain <select id="chain"></select></label>
<label>Pose <select id="pose"></select></label><label><input id="water" type="checkbox">Water</label>
<button id="reset">Reset view</button></header><div id="viewer"></div><div id="error" role="alert"></div>
<footer>Coordinates retained. A displayed pose is not a validated binding prediction.</footer>
<script src="3Dmol-min.js"></script><script id="data" type="application/json">__DATA__</script>
<script>
try {
 const data=JSON.parse(document.getElementById('data').textContent);
 const viewer=$3Dmol.createViewer('viewer',{backgroundColor:'white'});
 const receptor=viewer.addModel(data.structure,'cif');
 const poses=data.poses.map(p=>viewer.addModel(p,'sdf'));
 const chain=document.getElementById('chain'), pose=document.getElementById('pose');
 function option(select,value,label){const o=document.createElement('option');o.value=value;o.textContent=label;select.appendChild(o);}
 option(chain,'all','All chains');data.chains.forEach((c,i)=>option(chain,String(i),c||'(blank)'));
 option(pose,'none','Hide ligands');poses.forEach((p,i)=>option(pose,String(i),'Pose '+(i+1)));if(poses.length)pose.value='0';
 function draw(){
   receptor.setStyle({},{});
   const selection=chain.value==='all'?{}:{chain:data.chains[Number(chain.value)]};
   receptor.setStyle(selection,{cartoon:{colorscheme:'chain'},stick:{radius:0.12}});
   const water={...selection,resn:['HOH','WAT','DOD','H2O','SOL','TIP','TIP3']};
   receptor.setStyle(water,document.getElementById('water').checked?{sphere:{radius:0.3,color:'red'}}:{});
   poses.forEach((m,i)=>{m.setStyle({},{stick:{radius:0.22}});if(pose.value===String(i))m.show();else m.hide();});
   viewer.render();
 }
 if(data.box){const b=data.box;viewer.addBox({center:{x:b.center[0],y:b.center[1],z:b.center[2]},dimensions:{w:b.size[0],h:b.size[1],d:b.size[2]},wireframe:true,color:'#aa7733'});}
 [chain,pose,document.getElementById('water')].forEach(e=>e.addEventListener('change',draw));
 document.getElementById('reset').addEventListener('click',()=>{viewer.zoomTo();viewer.render();});
 draw();viewer.zoomTo();viewer.render();window.addEventListener('resize',()=>{viewer.resize();viewer.render();});
} catch(error){document.getElementById('error').textContent='Viewer could not start: '+error.message;}
</script></html>'''


def export_view(structure_file: Path, output_dir: Path, *, poses_file=None, site_file=None, model=None) -> dict:
    from rdkit import Chem
    output_dir = Path(output_dir)
    selection = select_structure(structure_file, output_dir / "structure.cif", model=model)
    data = {"structure": (output_dir / "structure.cif").read_text(), "poses": [], "box": None,
            "chains": [c.chain_id for c in read_structure(output_dir / "structure.cif").first_model.chains]}
    if poses_file:
        molecules = load_ligands(Path(poses_file))
        if len(molecules) > 200:
            raise ValueError("visualization export is limited to 200 explicitly selected poses")
        if any(not mol.GetNumConformers() for mol in molecules):
            raise ValueError("3D viewing requires existing ligand coordinates")
        with Chem.SDWriter(str(output_dir / "poses.sdf")) as writer:
            for mol in molecules:
                writer.write(mol)
                data["poses"].append(Chem.MolToMolBlock(mol))
    if site_file:
        data["box"] = BindingBox.from_dict(read_json(site_file)).to_dict()
    static = Path(__file__).resolve().parents[2] / "dashboard" / "static"
    for filename in ("3Dmol-min.js", "3Dmol-min.js.LICENSE.txt"):
        shutil.copyfile(static / filename, output_dir / filename)
    payload = json.dumps(data, allow_nan=False).replace("<", "\\u003c").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
    (output_dir / "view.html").write_text(_PAGE.replace("__DATA__", payload), encoding="utf-8")
    pymol = "from pathlib import Path\nfrom pymol import cmd\nroot = Path(__file__).resolve().parent\ncmd.load(str(root / 'structure.cif'), 'receptor')\ncmd.hide('everything', 'receptor')\ncmd.show('cartoon', 'receptor')\ncmd.show('sticks', 'receptor and hetatm and not solvent')\n"
    chimerax = "from pathlib import Path\nimport json\nfrom chimerax.core.commands import run\nroot = Path(__file__).resolve().parent\nrun(session, 'open ' + json.dumps(str(root / 'structure.cif')))\nrun(session, 'cartoon')\n"
    if data["poses"]:
        pymol += "cmd.load(str(root / 'poses.sdf'), 'poses')\ncmd.show('sticks', 'poses')\n"
        chimerax += "run(session, 'open ' + json.dumps(str(root / 'poses.sdf')))\n"
    pymol += "cmd.zoom()\n"
    (output_dir / "view_pymol.py").write_text(pymol)
    (output_dir / "view_chimerax.py").write_text(chimerax)
    report = {"schema_version": 1, "viewers": ["3Dmol", "PyMOL", "ChimeraX"],
              "pose_count": len(data["poses"]), "selection": selection,
              "external_viewers_launched": False, "coordinate_frame": "input; no alignment",
              "files": {p.name: file_hash(p) for p in output_dir.iterdir() if p.is_file()}}
    write_json(output_dir / "view_manifest.json", report)
    return report
