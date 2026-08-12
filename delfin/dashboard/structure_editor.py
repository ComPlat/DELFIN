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

The numbers are a layer over the molecule and nothing else. They read the
atoms the viewer already holds, so they travel with the cores through a drag,
an optimisation or a dynamic run; switching them off removes sprites and draws
a frame, and leaves the model -- and with it the bonds -- untouched.
"""

from __future__ import annotations


#: On-screen size of the numbers: the factor the high-resolution label texture
#: is down-scaled by, so a larger number stays sharp instead of blurring.
LABEL_SCALE_DEFAULT = 0.50

#: What a size control offers, smallest first. The smallest is what used to be
#: the largest but one: read at arm's length from a laptop screen, the old
#: ladder started far below useful and only its top rung could be read at all.
LABEL_SIZES = (('S', 0.28), ('M', 0.38), ('L', 0.50), ('XL', 0.66),
               ('XXL', 0.86))


def _atom_numbers_js():
    """The layer itself: ``window.__delfinAtomNumbers``.

    ``set(viewer, on, scale)`` switches the numbers on or off, ``refresh``
    brings them back onto the atoms after those have moved, ``setScale``
    resizes what is already there.
    """
    return (
        "window.__delfinAtomNumbers=(function(){\n"
        "  var RAD=0.5, RAD2=RAD*RAD, EPS=0.05, DEPTH_SIGN=1;\n"
        "  var DEFAULT_SCALE=%.3f;\n"
        # Every model in the viewer is numbered, each from its own zero: the
        # ORCA Builder's overlay shows two structures at once and both want
        # the numbering of the file they came from.
        "  function modelsOf(v){\n"
        "    var out=[];\n"
        "    try{\n"
        "      var ms=v.models;\n"
        "      if(ms&&ms.length){for(var i=0;i<ms.length;i++){if(ms[i])out.push(ms[i]);}}\n"
        "      else{var m=v.getModel();if(m)out.push(m);}\n"
        "    }catch(e){}\n"
        "    return out;\n"
        "  }\n"
        "  function atomCount(v){\n"
        "    var ms=modelsOf(v),n=0;\n"
        "    for(var i=0;i<ms.length;i++){\n"
        "      try{n+=(ms[i].selectedAtoms({})||[]).length;}catch(e){}\n"
        "    }\n"
        "    return n;\n"
        "  }\n"
        # Taking the sprites out by hand rather than with removeLabel: that
        # one draws a frame per label, so a hundred numbers meant a hundred
        # renders. Nothing here reaches the model, which is the point -- the
        # numbers must not cost the user their bonds.
        "  function clear(v){\n"
        "    if(!v)return;\n"
        "    var L=v.__delfinLbls||[];\n"
        "    for(var i=0;i<L.length;i++){\n"
        "      var lab=L[i]&&L[i].l;if(!lab)continue;\n"
        "      try{if(lab.sprite&&v.modelGroup)v.modelGroup.remove(lab.sprite);}catch(e){}\n"
        "      try{var k=(v.labels||[]).indexOf(lab);if(k>=0)v.labels.splice(k,1);}catch(e){}\n"
        "      try{lab.dispose();}catch(e){}\n"
        "    }\n"
        "    v.__delfinLbls=[];v.__delfinProj=[];\n"
        "  }\n"
        "  function hook(v){\n"
        "    if(v.__delfinLabelHooked)return;\n"
        "    v.__delfinLabelHooked=true;\n"
        "    var settle=function(){\n"
        "      try{if(refresh(v)&&typeof v.render==='function')v.render();}catch(e){}\n"
        "    };\n"
        "    var hs=v.__delfinInteractionEndHandlers||(v.__delfinInteractionEndHandlers=[]);\n"
        "    hs.push(settle);\n"
        "    var before=v.__delfinCleanup;\n"
        "    v.__delfinCleanup=function(){\n"
        "      var pos=hs.indexOf(settle);if(pos>=0)hs.splice(pos,1);\n"
        "      v.__delfinLbls=[];v.__delfinProj=[];\n"
        "      if(typeof before==='function'){try{before();}catch(e){}}\n"
        "    };\n"
        "  }\n"
        "  function build(v,scale){\n"
        "    if(!v||typeof v.addLabel!=='function')return 0;\n"
        "    clear(v);\n"
        "    if(scale!=null&&isFinite(+scale))v.__delfinLabelScale=+scale;\n"
        # alignment:center anchors the text box on its centre, so the number
        # stays on the atom centre at every zoom (default corner-anchoring
        # drifts aside as atoms shrink). Fall back to the string form if the
        # enum is unavailable.
        "    var al=(window.$3Dmol&&$3Dmol.SpriteAlignment&&$3Dmol.SpriteAlignment.center)\n"
        "      ?$3Dmol.SpriteAlignment.center:'center';\n"
        "    var ms=modelsOf(v),L=[];\n"
        "    for(var mi=0;mi<ms.length;mi++){\n"
        "      var atoms=[];try{atoms=ms[mi].selectedAtoms({})||[];}catch(e){atoms=[];}\n"
        "      for(var i=0;i<atoms.length;i++){\n"
        # fontSize 48 is a HIGH-RES texture kept sharp; the sprite is then
        # down-scaled in refresh() so the number appears small and crisp.
        # The fourth argument tells 3Dmol not to draw a frame per label.
        "        var a=atoms[i],lab=null;\n"
        "        try{lab=v.addLabel(String(i),{position:{x:a.x,y:a.y,z:a.z},\n"
        "          fontSize:48,fontColor:'black',alignment:al,\n"
        "          showBackground:false,inFront:true},undefined,true);}catch(e){lab=null;}\n"
        "        if(lab)L.push({a:a,l:lab});\n"
        "      }\n"
        "    }\n"
        "    v.__delfinLbls=L;v.__delfinProj=[];\n"
        "    if(L.length)hook(v);\n"
        "    return L.length;\n"
        "  }\n"
        # Called on every frame the molecule is drawn in, so a number is where
        # its atom is at that moment -- during a drag, an optimisation or a
        # dynamic run as much as after one. The atom object is held, not a
        # copy of its coordinates, which is why following it costs nothing.
        "  function refresh(v){\n"
        "    if(!v||v.__delfinDisposed)return false;\n"
        "    var L=v.__delfinLbls||[],m=L.length;\n"
        "    if(!m)return false;\n"
        "    if(atomCount(v)!==m&&!v.__delfinLabelRebuilding){\n"
        "      v.__delfinLabelRebuilding=true;\n"
        "      try{m=build(v,v.__delfinLabelScale);}finally{v.__delfinLabelRebuilding=false;}\n"
        "      L=v.__delfinLbls||[];\n"
        "      if(!m)return false;\n"
        "    }\n"
        "    var w;try{w=v.getView();}catch(e){return false;}\n"
        "    if(!w||w.length<8)return false;\n"
        "    var x=w[4],y=w[5],z=w[6],q=w[7];\n"
        "    var r11=1-2*(y*y+z*z), r12=2*(x*y-q*z), r13=2*(x*z+q*y);\n"
        "    var r21=2*(x*y+q*z), r22=1-2*(x*x+z*z), r23=2*(y*z-q*x);\n"
        "    var r31=2*(x*z-q*y), r32=2*(y*z+q*x), r33=1-2*(x*x+y*y);\n"
        "    var P=v.__delfinProj||(v.__delfinProj=[]);\n"
        "    var grid=Object.create(null);\n"
        "    for(var i=0;i<m;i++){\n"
        "      var a=L[i].a,sp=L[i].l&&L[i].l.sprite;\n"
        "      if(sp&&sp.position&&typeof sp.position.set==='function')\n"
        "        sp.position.set(a.x,a.y,a.z);\n"
        "      var p=P[i]||(P[i]=[0,0,0]);\n"
        "      p[0]=r11*a.x+r12*a.y+r13*a.z;\n"
        "      p[1]=r21*a.x+r22*a.y+r23*a.z;\n"
        "      p[2]=DEPTH_SIGN*(r31*a.x+r32*a.y+r33*a.z);\n"
        "      var gx=Math.floor(p[0]/RAD),gy=Math.floor(p[1]/RAD);\n"
        "      var key=gx+':'+gy,bucket=grid[key]||(grid[key]=[]);bucket.push(i);\n"
        "    }\n"
        # A number sits on top of its own atom (inFront) so it is never hidden
        # by it; one that is behind another atom is hidden here instead. The
        # grid keeps that to nearby candidates rather than every pair.
        "    for(var b0=0;b0<m;b0++){ var occ=false, pa=P[b0];\n"
        "      var agx=Math.floor(pa[0]/RAD),agy=Math.floor(pa[1]/RAD);\n"
        "      for(var ox=-1;ox<=1&&!occ;ox++){for(var oy=-1;oy<=1&&!occ;oy++){\n"
        "        var near=grid[(agx+ox)+':'+(agy+oy)]||[];\n"
        "        for(var n=0;n<near.length;n++){var b=near[n];if(b===b0)continue;var pb=P[b];\n"
        "          if(pb[2]>pa[2]+EPS){var ex=pa[0]-pb[0],ey=pa[1]-pb[1];\n"
        "            if(ex*ex+ey*ey<RAD2){occ=true;break;}}\n"
        "        }\n"
        "      }}\n"
        "      var s=L[b0].l&&L[b0].l.sprite; if(s) s.visible=!occ;\n"
        "    }\n"
        # The zoom factor f measures pixels-per-Angstrom in the SCREEN PLANE
        # at the rotation centre (world (0,0,rz), depth invariant under
        # rotation) along the screen-right model direction, so f tracks the
        # mouse-wheel zoom only, NOT rotation. DISP is the on-screen size the
        # toolbar's size control asked for.
        "    var DISP=(+v.__delfinLabelScale)||DEFAULT_SCALE, f=1;\n"
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
        "  function draw(v){\n"
        "    if(typeof v.render==='function'){try{v.render();}catch(e){}}\n"
        "  }\n"
        "  function setScale(v,scale){\n"
        "    if(!v||!isFinite(+scale))return false;\n"
        "    v.__delfinLabelScale=+scale;\n"
        "    if(!(v.__delfinLbls||[]).length)return false;\n"
        "    if(refresh(v))draw(v);\n"
        "    return true;\n"
        "  }\n"
        "  function set(v,on,scale){\n"
        "    if(!v)return 0;\n"
        "    if(!on){\n"
        "      var had=(v.__delfinLbls||[]).length;\n"
        "      clear(v);if(had)draw(v);\n"
        "      return 0;\n"
        "    }\n"
        "    var n=build(v,scale);\n"
        "    if(n){refresh(v);draw(v);}\n"
        "    return n;\n"
        "  }\n"
        "  return {set:set,build:build,clear:clear,refresh:refresh,setScale:setScale};\n"
        "})();"
    ) % LABEL_SCALE_DEFAULT


def atom_numbers_js():
    """The layer, installed once per page and only if it is not there yet."""
    return 'if(!window.__delfinAtomNumbers){\n' + _atom_numbers_js() + '\n}'


def show_atom_numbers_js(var='viewer', on=True, scale=None):
    """Number the atoms of the viewer held in the JS variable *var*.

    The numbers are read off the model the viewer already has, so this is the
    same call whether a molecule was just rendered or has been on screen for
    an hour, and it never needs the coordinates handed to it a second time.
    """
    size = LABEL_SCALE_DEFAULT if scale is None else float(scale)
    return (
        atom_numbers_js()
        + '\nwindow.__delfinAtomNumbers.set(%s,%s,%.3f);'
        % (var, 'true' if on else 'false', size)
    )
