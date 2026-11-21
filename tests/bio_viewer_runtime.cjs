const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

const source = fs.readFileSync(path.join(__dirname, '../delfin/dashboard/tab_agent.py'), 'utf8');
function extract(name) {
  const start = source.indexOf('function ' + name + '(');
  assert.ok(start >= 0, name);
  let depth = 0;
  const opening = source.indexOf('{', start);
  for (let index = opening; index < source.length; index++) {
    if (source[index] === '{') depth++;
    if (source[index] === '}' && --depth === 0) return source.slice(start, index + 1);
  }
  throw Error('Unclosed function ' + name);
}
const context = {};
vm.runInNewContext(['applyBioStyle', 'wireBio', 'hashOf', 'keyOf'].map(extract).join('\n'), context);
function viewer() {
  return {
    styles: [], additions: [], boxes: [], renders: 0,
    setStyle(selection, style) { this.styles.push({selection, style}); },
    addStyle(selection, style) { this.additions.push({selection, style}); },
    removeAllShapes() { this.boxes = []; },
    addBox(box) { this.boxes.push(box); },
    render() { this.renders++; }
  };
}
const chain = {value: '__all__'};
const waters = {checked: false};
const primary = viewer();
const card = {
  __delfinMol3DBio: {
    chains: ['A', 'B'], water_residue_names: ['HOH'],
    binding_site: {kind: 'binding_box', center: [1, 2, 3], size: [10, 12, 14]}
  },
  __delfinViewer: primary,
  querySelector(selector) { return selector === '.delfin-bio-chain' ? chain : waters; }
};
context.wireBio(card, primary);
assert.equal(primary.styles.filter(entry => entry.style.cartoon).length, 2);
assert.equal(primary.styles.filter(entry => entry.style.sphere).length, 0);
assert.equal(primary.boxes[0].dimensions.w, 10);
chain.value = 'B';
waters.checked = true;
chain.onchange();
assert.equal(card.__delfinMol3DBioState.chain, 'B');
assert.equal(primary.styles.at(-1).selection.chain, 'B');
assert.ok(primary.styles.at(-1).style.sphere);
assert.equal(primary.boxes.length, 1);
const fullscreen = viewer();
const fullscreenChain = {value: ''};
const fullscreenWaters = {checked: false};
context.wireBio(card, fullscreen, {
  querySelector(selector) { return selector === '.delfin-bio-chain' ? fullscreenChain : fullscreenWaters; }
});
assert.equal(fullscreenChain.value, 'B');
assert.equal(fullscreenWaters.checked, true);
fullscreenChain.value = 'A';
fullscreenChain.onchange();
assert.equal(chain.value, 'A');
assert.ok(primary.renders > 0);
card.__delfinMol3DBio.binding_site = {kind: 'residue_selection'};
card.__delfinMol3DBio.binding_site_serials = [42, 43];
context.applyBioStyle(primary, card.__delfinMol3DBio, {chain: 'A', waters: false});
assert.equal(primary.additions.at(-1).selection.serial[0], 42);
assert.equal(primary.additions.at(-1).selection.chain, 'A');
const attributes = {'data-mol3d': 'same coordinates', 'data-mol3d-bio': 'box one', 'data-mol3d-fmt': 'pdb'};
const cacheCard = {getAttribute(name) { return attributes[name]; }};
const firstKey = context.keyOf({id: 'same-stage'}, cacheCard);
attributes['data-mol3d-bio'] = 'box two';
assert.notEqual(context.keyOf({id: 'same-stage'}, cacheCard), firstKey);
console.log('PASS: cartoon chains, hidden/toggled waters, binding boxes, exact atom selections, fullscreen state and cache invalidation.');
