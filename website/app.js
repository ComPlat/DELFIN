/* Buildless, dependency-free interactions. Content is in content.js. */
(() => {
  'use strict';
  const content = window.DELFIN_CONTENT;
  if (!content) return;
  document.documentElement.classList.add('js');
  const $ = (selector) => document.querySelector(selector);
  const reducedMotion = matchMedia('(prefers-reduced-motion: reduce)');
  let paused = reducedMotion.matches;
  let manualMotionChoice = false;
  document.querySelectorAll('[data-content]').forEach((node) => {
    const value = node.dataset.content.split('.').reduce((value, key) => value?.[key], content);
    if (typeof value === 'string') node.textContent = value;
  });

  function animatePanel(node) {
    node.classList.remove('panel-enter');
    if (paused) return;
    void node.offsetWidth;
    node.classList.add('panel-enter');
  }
  function choiceButtons(items, parent, className, render, onSelect) {
    items.forEach((item, index) => {
      const button = document.createElement('button');
      button.type = 'button';
      button.className = className;
      button.innerHTML = render(index);
      const text = document.createTextNode(item.label || item.name);
      button.insertBefore(text, button.lastElementChild?.tagName === 'B' ? button.lastElementChild : null);
      button.setAttribute('aria-pressed', String(index === 0));
      button.addEventListener('click', () => {
        parent.querySelectorAll('button').forEach((node) => node.setAttribute('aria-pressed', String(node === button)));
        onSelect(index);
      });
      parent.append(button);
    });
    onSelect(0);
  }
  choiceButtons(content.workflow, $('.workflow-steps'), 'workflow-step', (i) => `<span>0${i + 1}</span>`, (index) => {
    const item = content.workflow[index];
    for (const key of ['tag', 'title', 'text', 'detail', 'code']) $('#flow-' + key).textContent = item[key];
    animatePanel($('.workflow-body'));
  });
  choiceButtons(content.research, $('.research-choices'), 'research-choice', (i) => `<span>0${i + 1}</span><b aria-hidden="true">↗</b>`, (index) => {
    const item = content.research[index];
    $('#research-symbol').textContent = item.symbol;
    $('#research-note').textContent = item.note;
    $('#research-heading').textContent = item.title;
    $('#research-text').textContent = item.text;
    $('#research-tags').replaceChildren(...item.tags.map((tag) => {
      const node = document.createElement('span'); node.textContent = tag; return node;
    }));
    animatePanel($('.research-panel'));
  });
  content.animals.forEach((animal) => {
    const card = document.createElement('a');
    card.className = 'animal-card reveal';
    card.href = animal.href;
    const img = document.createElement('img');
    img.src = animal.image; img.alt = animal.name + ' logo'; img.width = 100; img.height = 100; img.loading = 'lazy';
    const title = document.createElement('h3'); title.textContent = animal.name;
    const role = document.createElement('p'); role.className = 'animal-role'; role.textContent = animal.role;
    const description = document.createElement('p'); description.textContent = animal.description;
    const arrow = document.createElement('span'); arrow.className = 'animal-arrow'; arrow.textContent = '↗'; arrow.setAttribute('aria-hidden', 'true');
    card.append(img, title, role, description, arrow);
    $('.animal-grid').append(card);
  });

  // Real videos only; nothing is fetched from a video provider until clicked.
  const dialog = $('#video-dialog');
  let videoTrigger;
  (content.videos || []).forEach((video) => {
    $('#watch').hidden = false;
    const button = document.createElement('button');
    button.className = 'video-card reveal'; button.type = 'button';
    const poster = document.createElement('img'); poster.src = video.poster; poster.alt = ''; poster.loading = 'lazy';
    const title = document.createElement('span'); title.textContent = '▶  ' + video.title;
    const description = document.createElement('p'); description.textContent = video.description || '';
    button.append(poster, title, description);
    button.addEventListener('click', () => {
      videoTrigger = button;
      $('#video-dialog-title').textContent = video.title;
      let player;
      if (video.type === 'youtube') {
        const url = new URL(video.src);
        if (url.protocol !== 'https:' || url.hostname !== 'www.youtube-nocookie.com' || !url.pathname.startsWith('/embed/')) return;
        player = document.createElement('iframe');
        player.src = url.href; player.title = video.title;
        player.allow = 'fullscreen; picture-in-picture'; player.allowFullscreen = true;
        player.referrerPolicy = 'strict-origin-when-cross-origin';
      } else {
        player = document.createElement('video'); player.src = video.src; player.controls = true; player.playsInline = true;
        player.poster = video.poster; player.preload = 'metadata';
        if (video.captions) {
          const track = document.createElement('track'); track.kind = 'captions'; track.src = video.captions;
          track.srclang = video.language || 'en'; track.label = video.captionLabel || 'English'; track.default = true; player.append(track);
        }
      }
      $('#video-player').replaceChildren(player); dialog.showModal(); document.body.classList.add('dialog-open');
    });
    $('.video-grid').append(button);
  });
  $('#close-video').addEventListener('click', () => dialog.close());
  dialog.addEventListener('click', (event) => { if (event.target === dialog) dialog.close(); });
  dialog.addEventListener('close', () => {
    $('#video-player').replaceChildren(); document.body.classList.remove('dialog-open'); videoTrigger?.focus();
  });

  $('#copy-install').addEventListener('click', async () => {
    try {
      await navigator.clipboard.writeText('pip install delfin-complat');
      $('#copy-status').textContent = 'Copied. Ready for your terminal.';
    } catch {
      const range = document.createRange(); range.selectNodeContents($('.install-command code'));
      const selection = window.getSelection(); selection.removeAllRanges(); selection.addRange(range);
      $('#copy-status').textContent = 'Command selected. Press Ctrl+C or ⌘C to copy.';
    }
  });

  const revealObserver = 'IntersectionObserver' in window ? new IntersectionObserver((entries) => {
    entries.forEach((entry) => { if (entry.isIntersecting) { entry.target.classList.add('visible'); revealObserver.unobserve(entry.target); } });
  }, {threshold: 0.08}) : null;
  document.querySelectorAll('.reveal').forEach((node) => {
    if (revealObserver) revealObserver.observe(node); else node.classList.add('visible');
  });
  let scrollPending = false;
  function updateProgress() {
    const height = document.documentElement.scrollHeight - innerHeight;
    $('.reading-progress').style.transform = `scaleX(${height > 0 ? scrollY / height : 0})`;
    scrollPending = false;
  }
  addEventListener('scroll', () => { if (!scrollPending) { scrollPending = true; requestAnimationFrame(updateProgress); } }, {passive:true});
  updateProgress();

  // Actual MANTA XYZ coordinates and inferred display connectivity, in molecule.js.
  // Canvas rendering also works on cluster browsers without WebGL acceleration.
  const canvas = $('#molecule');
  const ctx = canvas.getContext('2d');
  const model = window.DELFIN_MOLECULE;
  if (!model || !ctx) {
    $('.structure-source').textContent = 'The structure viewer is unavailable. Download the XYZ coordinates to inspect the molecule.';
    return;
  }
  let width=1,height=1,last=0,frame=null,visible=true;
  let rotationY=.55,rotationX=-.3,zoom=1,showHydrogens=true,dragging=false,lastPointer=null;
  const atoms=model.atoms;
  const bonds=model.bonds;
  const center=atoms.reduce((s,a)=>({x:s.x+a.x/atoms.length,y:s.y+a.y/atoms.length,z:s.z+a.z/atoms.length}),{x:0,y:0,z:0});
  const centered=atoms.map(a=>({...a,x:a.x-center.x,y:a.y-center.y,z:a.z-center.z}));
  const extent=Math.max(...centered.map(a=>Math.hypot(a.x,a.y,a.z)));
  const palette={Ir:['#fff1ca','#cfaa6a','#80602f'],N:['#bfdcff','#5684cc','#304c7a'],C:['#c8d4df','#7e92a5','#42576b'],H:['#fff','#e4eaf0','#a4b3c1']};
  function project(atom){
    const cy=Math.cos(rotationY),sy=Math.sin(rotationY),cx=Math.cos(rotationX),sx=Math.sin(rotationX);
    const x=atom.x*cy+atom.z*sy,z=-atom.x*sy+atom.z*cy;
    const y=atom.y*cx-z*sx,depth=atom.y*sx+z*cx;
    const scale=Math.min(width,height)*.41/extent*zoom,perspective=40/(40+depth);
    return {...atom,x:width/2+x*scale*perspective,y:height/2+y*scale*perspective,z:depth,scale:scale*perspective,r:(atom.element==='Ir'?.48:atom.element==='H'?.14:.27)*scale*perspective};
  }
  function draw(){
    ctx.clearRect(0,0,width,height);
    const projected=centered.map(project);
    const objects=[...projected.filter(a=>showHydrogens||a.element!=='H').map(atom=>({type:'atom',z:atom.z,atom})),...bonds.filter(([a,b])=>showHydrogens||(atoms[a].element!=='H'&&atoms[b].element!=='H')).map(([a,b])=>({type:'bond',z:(projected[a].z+projected[b].z)/2,a:projected[a],b:projected[b]}))].sort((a,b)=>b.z-a.z);
    objects.forEach(object=>{
      if(object.type==='bond'){
        const {a,b}=object;
        const gradient=ctx.createLinearGradient(a.x,a.y,b.x,b.y);
        gradient.addColorStop(0,(palette[a.element]||palette.C)[1]);gradient.addColorStop(1,(palette[b.element]||palette.C)[1]);
        ctx.lineCap='round';ctx.strokeStyle=gradient;ctx.lineWidth=Math.max(1.5,(a.scale+b.scale)*.075);
        ctx.beginPath();ctx.moveTo(a.x,a.y);ctx.lineTo(b.x,b.y);ctx.stroke();
        ctx.strokeStyle='#ffffff25';ctx.lineWidth=Math.max(.5,(a.scale+b.scale)*.018);ctx.stroke();
      }else{
        const {x,y,r,element}=object.atom;
        const gradient=ctx.createRadialGradient(x-r*.35,y-r*.4,r*.05,x,y,r);
        const colors=palette[element]||palette.C;
        gradient.addColorStop(0,colors[0]);gradient.addColorStop(.55,colors[1]);gradient.addColorStop(1,colors[2]);
        ctx.fillStyle=gradient;ctx.beginPath();ctx.arc(x,y,r,0,Math.PI*2);ctx.fill();
      }
    });
  }
  function tick(now){
    frame=null;
    if(paused||!visible||document.hidden){last=0;return;}
    if(last&&!dragging)rotationY+=Math.min((now-last)/1000,.05)*.12;
    last=now;draw();frame=requestAnimationFrame(tick);
  }
  function startAnimation(){if(frame===null&&!paused&&visible&&!document.hidden)frame=requestAnimationFrame(tick);}
  function resize(){
    const box=canvas.getBoundingClientRect();width=box.width;height=box.height;
    const dpr=Math.min(devicePixelRatio||1,2);canvas.width=Math.round(width*dpr);canvas.height=Math.round(height*dpr);
    ctx.setTransform(dpr,0,0,dpr,0,0);draw();
  }
  addEventListener('resize',resize,{passive:true});resize();
  canvas.addEventListener('pointerdown',event=>{
    if(event.button!==0)return;
    dragging=true;lastPointer={x:event.clientX,y:event.clientY};canvas.setPointerCapture(event.pointerId);
  });
  canvas.addEventListener('pointermove',event=>{
    if(!dragging)return;
    rotationY+=(event.clientX-lastPointer.x)*.009;rotationX+=(event.clientY-lastPointer.y)*.009;
    lastPointer={x:event.clientX,y:event.clientY};draw();
  });
  const release=()=>{dragging=false;lastPointer=null;};
  canvas.addEventListener('pointerup',release);canvas.addEventListener('pointercancel',release);canvas.addEventListener('lostpointercapture',release);
  canvas.addEventListener('keydown',event=>{
    if(!['ArrowLeft','ArrowRight','ArrowUp','ArrowDown','+','=','-','0'].includes(event.key))return;
    event.preventDefault();
    if(event.key==='ArrowLeft')rotationY-=.15;if(event.key==='ArrowRight')rotationY+=.15;
    if(event.key==='ArrowUp')rotationX-=.15;if(event.key==='ArrowDown')rotationX+=.15;
    if(event.key==='+'||event.key==='=')zoom=Math.min(1.5,zoom+.1);if(event.key==='-')zoom=Math.max(.65,zoom-.1);
    if(event.key==='0'){rotationY=.55;rotationX=-.3;zoom=1;}draw();
  });
  $('#molecule-reset').addEventListener('click',()=>{rotationY=.55;rotationX=-.3;zoom=1;draw();});
  $('#molecule-hydrogens').addEventListener('click',()=>{
    showHydrogens=!showHydrogens;$('#molecule-hydrogens').setAttribute('aria-pressed',String(showHydrogens));draw();
  });
  if('IntersectionObserver' in window)new IntersectionObserver(([entry])=>{visible=entry.isIntersecting;if(visible)startAnimation();},{threshold:0}).observe(canvas);
  document.addEventListener('visibilitychange',()=>{if(!document.hidden)startAnimation();});
  function setMotion(){
    document.body.classList.toggle('motion-paused',paused);
    $('.motion-toggle').textContent=paused?'Enable motion':'Pause motion';
    $('.motion-toggle').setAttribute('aria-pressed',String(paused));
    if(paused){if(frame!==null)cancelAnimationFrame(frame);frame=null;last=0;draw();}else startAnimation();
  }
  $('.motion-toggle').addEventListener('click',()=>{manualMotionChoice=true;paused=!paused;setMotion();});
  reducedMotion.addEventListener('change',()=>{if(!manualMotionChoice){paused=reducedMotion.matches;setMotion();}});
  setMotion();
})();
