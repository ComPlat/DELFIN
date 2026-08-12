"""The structure editor, as a part rather than as a tab.

Everything the Submit tab's molecule preview can do -- picking, manipulating,
drawing, the browser force fields and the ones that run on the server, the
held internal coordinates, the polyhedra, the numbers on the atoms -- belongs
to a viewer, not to one tab. This module is where it lives so a second tab can
have the same editor over its own coordinates.

It is being assembled a piece at a time, and each piece has to leave the Submit
tab behaving exactly as it did. The first is the atom numbering, which was
reachable only from the ORCA Builder: its size control resolved the viewer as
``window._orcaBuildViewer``, one global for one tab. Numbering is a property of
a viewer, so it takes the viewer it is meant for and any tab can show it.
"""

from __future__ import annotations


#: On-screen size of the numbers: the factor the high-resolution label texture
#: is down-scaled by, so a larger number stays sharp instead of blurring.
LABEL_SCALE_DEFAULT = 0.20

#: What a size control offers, smallest first.
LABEL_SIZES = (('S', 0.11), ('M', 0.15), ('L', 0.20), ('XL', 0.28),
               ('XXL', 0.38))


_LABEL_DEPTH_PATCH_JS = (
    "(function(){\n"
    "  var v=__VAR__, RAD=0.5, RAD2=RAD*RAD, EPS=0.05, DEPTH_SIGN=1;\n"
    "  var disposed=false, updateTimer=0;\n"
    "  if(v.__delfinLabelDepthPatched) return;\n"
    "  v.__delfinLabelDepthPatched=true;\n"
    "  function update(){\n"
    "    if(disposed||v.__delfinDisposed||v.__delfinInteracting) return false;\n"
    "    var w; try{w=v.getView();}catch(e){return;}\n"
    "    if(!w||w.length<8) return false;\n"
    "    var x=w[4],y=w[5],z=w[6],q=w[7];\n"
    "    var r11=1-2*(y*y+z*z), r12=2*(x*y-q*z), r13=2*(x*z+q*y);\n"
    "    var r21=2*(x*y+q*z), r22=1-2*(x*x+z*z), r23=2*(y*z-q*x);\n"
    "    var r31=2*(x*z-q*y), r32=2*(y*z+q*x), r33=1-2*(x*x+y*y);\n"
    "    var L=v.__delfinLbls||[], m=L.length;\n"
    "    var P=v.__delfinProj||(v.__delfinProj=[]);\n"
    "    var grid=Object.create(null);\n"
    "    for(var i=0;i<m;i++){ var c=L[i].c; var p=P[i]||(P[i]=[0,0,0]);\n"
    "      p[0]=r11*c[0]+r12*c[1]+r13*c[2];\n"
    "      p[1]=r21*c[0]+r22*c[1]+r23*c[2];\n"
    "      p[2]=DEPTH_SIGN*(r31*c[0]+r32*c[1]+r33*c[2]);\n"
    "      var gx=Math.floor(p[0]/RAD),gy=Math.floor(p[1]/RAD);\n"
    "      var key=gx+':'+gy,bucket=grid[key]||(grid[key]=[]);bucket.push(i);\n"
    "    }\n"
    "    for(var a=0;a<m;a++){ var occ=false, pa=P[a];\n"
    "      var agx=Math.floor(pa[0]/RAD),agy=Math.floor(pa[1]/RAD);\n"
    "      for(var ox=-1;ox<=1&&!occ;ox++){for(var oy=-1;oy<=1&&!occ;oy++){\n"
    "        var near=grid[(agx+ox)+':'+(agy+oy)]||[];\n"
    "        for(var n=0;n<near.length;n++){var b=near[n];if(b===a)continue;var pb=P[b];\n"
    "          if(pb[2]>pa[2]+EPS){var ex=pa[0]-pb[0],ey=pa[1]-pb[1];\n"
    "            if(ex*ex+ey*ey<RAD2){occ=true;break;}}\n"
    "        }\n"
    "      }}\n"
    "      var s=L[a].l&&L[a].l.sprite; if(s) s.visible=!occ;\n"
    "    }\n"
    "    // Label sizing.  The text texture is rendered at a high fontSize\n"
    "    // (crisp), and we always DOWN-scale the sprite (DISP) so it looks\n"
    "    // small and sharp -- never up-scaled (which is what caused blur).\n"
    "    // The zoom factor f measures pixels-per-Angstrom in the SCREEN\n"
    "    // PLANE at the rotation centre (world (0,0,rz), depth invariant\n"
    "    // under rotation) along the screen-right model direction, so f\n"
    "    // tracks the mouse-wheel zoom only, NOT rotation.  DISP is the base\n"
    "    // on-screen size, set from the toolbar's size control.\n"
    "    var DISP=(+v.__delfinLabelScale)||__DISP__, f=1;\n"
    "    if(typeof v.modelToScreen==='function'&&m>0){ try{\n"
    "      var pcx=-w[0],pcy=-w[1],pcz=-w[2];\n"
    "      var q1=v.modelToScreen({x:pcx,y:pcy,z:pcz});\n"
    "      var q2=v.modelToScreen({x:pcx+r11,y:pcy+r12,z:pcz+r13});\n"
    "      var ppa=Math.sqrt((q1.x-q2.x)*(q1.x-q2.x)+(q1.y-q2.y)*(q1.y-q2.y));\n"
    "      if(ppa>0){ var base=v.__delfinPPA0||(v.__delfinPPA0=ppa);\n"
    "        f=ppa/base; if(f<0.4)f=0.4; if(f>2.5)f=2.5; }\n"
    "    }catch(e){} }\n"
    "    for(var g=0;g<m;g++){ var sg=L[g].l&&L[g].l.sprite;\n"
    "      if(!sg||!sg.scale) continue;\n"
    "      if(!L[g].s0&&sg.scale.x>0) L[g].s0=[sg.scale.x,sg.scale.y];\n"
    "      if(L[g].s0){ sg.scale.x=L[g].s0[0]*DISP*f; sg.scale.y=L[g].s0[1]*DISP*f; }\n"
    "    }\n"
    "    return true;\n"
    "  }\n"
    "  function applySettled(){\n"
    "    updateTimer=0;\n"
    "    if(disposed||v.__delfinDisposed) return;\n"
    "    try{if(update()&&typeof v.render==='function')v.render();}catch(e){}\n"
    "  }\n"
    "  function schedule(delay){\n"
    "    if(disposed||v.__delfinDisposed)return;\n"
    "    if(updateTimer)clearTimeout(updateTimer);\n"
    "    updateTimer=setTimeout(applySettled,delay||0);\n"
    "  }\n"
    "  var afterInteraction=function(){schedule(0);};\n"
    "  var handlers=v.__delfinInteractionEndHandlers||(v.__delfinInteractionEndHandlers=[]);\n"
    "  handlers.push(afterInteraction);\n"
    "  var previousCleanup=v.__delfinCleanup;\n"
    "  v.__delfinCleanup=function(){\n"
    "    if(disposed)return;disposed=true;\n"
    "    if(updateTimer)clearTimeout(updateTimer);updateTimer=0;\n"
    "    var hs=v.__delfinInteractionEndHandlers||[];\n"
    "    var pos=hs.indexOf(afterInteraction);if(pos>=0)hs.splice(pos,1);\n"
    "    if(typeof previousCleanup==='function'){try{previousCleanup();}catch(e){}}\n"
    "  };\n"
    "  schedule(0);\n"
    "})();"
)


def label_scale_setter_js():
    """The one-off script that teaches the page how to resize numbers.

    It used to be written inside the depth patch and resolved the viewer as
    ``window._orcaBuildViewer``, so only the ORCA Builder could ever call it.
    It takes the viewer now, and every tab passes its own.
    """
    return (
        "window.__delfinSetLabelScale=function(scale, viewer){"
        "if(!viewer) return false;"
        "viewer.__delfinLabelScale=scale;"
        "var hs=viewer.__delfinInteractionEndHandlers||[];"
        "for(var i=0;i<hs.length;i++){try{hs[i]();}catch(e){}}"
        "return true;};"
    )


def atom_number_labels_js(full_xyz, var='viewer', scale=None):
    """Return JS fragment adding atom-index labels with occlusion culling.

    Numbers sit at exact atom centres and on top (inFront:true) so they never
    drift with zoom and an atom never hides its own number. After camera
    interaction settles, _LABEL_DEPTH_PATCH_JS hides numbers of atoms that
    are behind other atoms without adding work to mouse-move renders."""
    size = LABEL_SCALE_DEFAULT if scale is None else float(scale)
    lines = full_xyz.split('\n')
    try:
        n_atoms = int(lines[0].strip())
    except (ValueError, IndexError):
        return ''
    # alignment:center anchors the text box on its centre, so the number
    # stays on the atom centre at every zoom (default corner-anchoring drifts
    # aside as atoms shrink).  Fall back to the string form if the enum is
    # unavailable.
    preamble = [
        'var __delfinAL=($3Dmol&&$3Dmol.SpriteAlignment&&$3Dmol.SpriteAlignment.center)'
        '?$3Dmol.SpriteAlignment.center:"center";',
        '%s.__delfinLbls=%s.__delfinLbls||[];' % (var, var),
        # carry the size chosen in the toolbar into every new viewer
        '%s.__delfinLabelScale=%.3f;' % (var, size),
    ]
    pushes = []
    for i, line in enumerate(lines[2: 2 + n_atoms]):
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            continue
        pushes.append(
            # fontSize 48 = HIGH-RES texture (kept sharp); the sprite is
            # then down-scaled (DISP) in _LABEL_DEPTH_PATCH_JS so the number
            # appears small and crisp.  alignment centred; no background.
            '%s.__delfinLbls.push({c:[%.6f,%.6f,%.6f],l:%s.addLabel("%d",'
            '{position:{x:%.6f,y:%.6f,z:%.6f},fontSize:48,fontColor:"black",'
            'alignment:__delfinAL,showBackground:false,inFront:true})});'
            % (var, x, y, z, var, i, x, y, z)
        )
    if not pushes:
        return ''
    depth_patch = (
        _LABEL_DEPTH_PATCH_JS
        .replace('__VAR__', var)
        .replace('__DISP__', '%.3f' % size)
    )
    return '\n    '.join(
        preamble + pushes + [depth_patch]
    )
