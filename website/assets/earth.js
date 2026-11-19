(() => {
  'use strict';
  const canvas = document.querySelector('.purpose-globe canvas');
  const context = canvas?.getContext('2d');
  if (!context) return;
  const motionPreference = matchMedia('(prefers-reduced-motion: reduce)');
  const map = new Image();
  map.addEventListener('load', () => {
    const texture = document.createElement('canvas');
    texture.width = 720;
    texture.height = 360;
    const textureContext = texture.getContext('2d');
    if (!textureContext) return;
    textureContext.drawImage(map, 0, 0, 720, 360);
    const pixels = textureContext.getImageData(0, 0, 720, 360).data;
    const size = 320;
    canvas.width = size;
    canvas.height = size;
    const surface = context.createImageData(size, size);
    const samples = [];
    const tilt = 12 * Math.PI / 180;
    for (let row = 0; row < size; row++) {
      for (let column = 0; column < size; column++) {
        const horizontal = (column + .5 - size / 2) / (size / 2);
        const vertical = (size / 2 - row - .5) / (size / 2);
        const distance = horizontal ** 2 + vertical ** 2;
        if (distance > 1) continue;
        const depth = Math.sqrt(1 - distance);
        const latitude = Math.asin(vertical * Math.cos(tilt) + depth * Math.sin(tilt));
        const longitude = Math.atan2(horizontal, depth * Math.cos(tilt) - vertical * Math.sin(tilt));
        const light = .65 + .35 * Math.max(0, -.4 * horizontal + .3 * vertical + .86 * depth);
        samples.push({offset: (row * size + column) * 4, latitude, longitude, light,
          alpha: Math.min(1, (1 - Math.sqrt(distance)) * size / 2) * 255});
      }
    }
    let rotation = 15;
    let frame = null;
    let last = 0;
    let visible = false;
    function draw() {
      samples.forEach(({offset, latitude, longitude, light, alpha}) => {
        const degrees = longitude * 180 / Math.PI + rotation;
        const textureColumn = Math.floor(((degrees + 180) % 360 + 360) % 360 * 2);
        const textureRow = Math.min(359, Math.floor((90 - latitude * 180 / Math.PI) * 2));
        const source = (textureRow * 720 + textureColumn) * 4;
        const latitudeLine = Math.abs((latitude * 180 / Math.PI + 90) % 30 - 15);
        const longitudeLine = Math.abs(((degrees % 30) + 30) % 30 - 15);
        const grid = latitudeLine > 14.65 || longitudeLine > 14.65 ? .13 : 0;
        for (let channel = 0; channel < 3; channel++) {
          surface.data[offset + channel] = pixels[source + channel] * light * (1 - grid) + 210 * grid;
        }
        surface.data[offset + 3] = alpha;
      });
      context.putImageData(surface, 0, 0);
    }
    function running() {
      return visible && !document.hidden && !document.body.classList.contains('motion-paused') &&
        (!motionPreference.matches || document.querySelector('.motion-toggle')?.getAttribute('aria-pressed') === 'false');
    }
    function tick(now) {
      frame = null;
      if (!running()) { last = 0; return; }
      if (!last || now - last >= 40) {
        if (last) rotation = (rotation + Math.min(now - last, 100) * .004) % 360;
        last = now;
        draw();
      }
      frame = requestAnimationFrame(tick);
    }
    function syncMotion() {
      if (running()) {
        if (frame === null) frame = requestAnimationFrame(tick);
      } else {
        if (frame !== null) cancelAnimationFrame(frame);
        frame = null;
        last = 0;
      }
    }
    draw();
    canvas.classList.add('ready');
    if ('IntersectionObserver' in window) {
      new IntersectionObserver(([entry]) => {
        visible = entry.isIntersecting;
        syncMotion();
      }).observe(canvas);
    } else {
      visible = true;
    }
    new MutationObserver(syncMotion).observe(document.body, {attributes: true, attributeFilter: ['class']});
    document.addEventListener('visibilitychange', syncMotion);
    motionPreference.addEventListener('change', syncMotion);
    syncMotion();
  });
  map.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(
    new XMLSerializer().serializeToString(document.querySelector('.purpose-globe svg'))
  );
})();
