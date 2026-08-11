"""One way of putting a large text file on screen, shared by every viewer.

A file used to arrive in the page as a single text node under
``white-space: pre-wrap``. Before the browser can paint anything it then has
to line-break and lay out every line of it. Measured in headless Chromium on
a 3.3 MB, 73 479 line ORCA output -- the shape of an ``S1.out`` -- that is
7.5 seconds during which the main thread answers nothing, which is what the
browser's "page unresponsive" prompt is reporting. Loading a smaller window
of the file does not help much: the cost is per line on screen, and a window
big enough to be worth scrolling is still tens of thousands of lines.

Here the same text is cut into blocks of whole lines, and every block
declares up front how tall it is going to be. ``content-visibility: auto``
then lets the browser skip layout and paint for every block that is off
screen, and it only ever does the work for the two or three blocks that are
actually in view. The same 3.3 MB takes about 20 ms.

Nothing about the text changes. The blocks concatenate back to exactly the
string that went in -- ``box.textContent === text`` -- so selection, the
browser's own Ctrl+F, and every character offset the search machinery works
in still address the whole chunk, not a window of it.

The block heights are predicted rather than measured, because measuring is
the expensive thing we are avoiding. For a monospace font that prediction is
exact: measure one line's height and one character's width, then count the
rows each line wraps into. Against the real layout of that same file the
predicted total came out to the pixel.
"""

import json

# Blocks end at whichever bound comes first, and never inside a line: a block
# boundary is a forced line break, so splitting a line across two blocks would
# change what the file looks like. The line bound keeps blocks cheap to
# re-render for a highlight; the character bound covers files whose lines are
# long enough that 400 of them would be a lot of layout in one block.
TEXT_VIEW_BLOCK_LINES = 400
TEXT_VIEW_BLOCK_CHARS = 24000


def text_view_css():
    """Styles for the block viewer. Safe to include more than once."""
    return (
        # The look of the old single text node, moved out of the inline style
        # so that the hidden probe which measures the font can wear the same
        # class and get the same answer.
        '.dtv-text { white-space: pre-wrap; overflow-wrap: anywhere;'
        ' word-break: break-word; font-family: monospace; font-size: 12px;'
        ' line-height: 1.3; }'
        # The whole point: a block that is off screen costs nothing.
        '.dtv-blk { content-visibility: auto; }'
        '.dtv-probe { position: absolute; visibility: hidden; height: auto;'
        ' top: 0; left: -99999px; pointer-events: none; }'
        '.dtv-mark { background: #fff59d; padding: 0 2px; }'
        '.dtv-mark.dtv-current { background: #ffcc80; }'
    )


def text_view_bootstrap_js():
    """The runtime, as one idempotent script.

    Every tab that shows text registers this, and the Calculations browser is
    built three times over (Calculations, Archive, Office), so it has to be
    fine to run repeatedly.
    """
    return _BOOTSTRAP_JS.replace(
        '__BLOCK_LINES__', str(TEXT_VIEW_BLOCK_LINES)
    ).replace(
        '__BLOCK_CHARS__', str(TEXT_VIEW_BLOCK_CHARS)
    )


def set_text_js(el_expr, text):
    """JS that fills the element named by ``el_expr`` with ``text``.

    Carries the runtime with it. A viewer can be asked to show a file before
    the page's startup scripts have run -- and if the runtime were missing at
    that moment the file would simply not appear -- so the definition travels
    with the first thing that needs it. It is a no-op once defined.
    """
    return (
        text_view_bootstrap_js()
        + f'window.__delfinTextView.setText({el_expr}, {json.dumps(text)});'
    )


_BOOTSTRAP_JS = r"""
(function () {
    if (window.__delfinTextView) return;

    var BLOCK_LINES = __BLOCK_LINES__;
    var BLOCK_CHARS = __BLOCK_CHARS__;

    function escapeHtml(s) {
        return s.replace(/[&<>"]/g, function (c) {
            return {'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'}[c];
        });
    }

    /* One line's height and one character's width, for the width the text is
       being shown at. Both come from a hidden copy wearing the same class, so
       a change to the stylesheet cannot make the prediction drift away from
       what is on screen. Cached until the column is resized. */
    function metrics(el) {
        var width = el.clientWidth || 0;
        var cached = el.__dtvMetrics;
        if (cached && cached.width === width && width > 0) return cached;

        var probe = document.createElement('div');
        probe.className = (el.className || '') + ' dtv-probe';
        probe.style.width = (width > 0 ? width : 600) + 'px';
        var rows = [];
        for (var i = 0; i < 11; i++) rows.push('x');
        probe.textContent = rows.join('\n');
        (el.parentElement || document.body).appendChild(probe);
        var lineHeight = probe.getBoundingClientRect().height / 11;

        /* An inline span in a wrapping block wraps too, and would measure the
           column rather than the glyph. Stop the wrapping for this one read. */
        probe.textContent = '';
        probe.style.whiteSpace = 'pre';
        var span = document.createElement('span');
        span.textContent = new Array(201).join('x');
        probe.appendChild(span);
        var charWidth = span.getBoundingClientRect().width / 200;
        probe.remove();

        if (!(lineHeight > 0)) lineHeight = 16;
        if (!(charWidth > 0)) charWidth = 7;
        var m = {
            width: width,
            lineHeight: lineHeight,
            charWidth: charWidth,
            cols: Math.max(1, Math.floor((width > 0 ? width : 600) / charWidth))
        };
        el.__dtvMetrics = m;
        return m;
    }

    /* How many rows a line wraps into: greedy at spaces, and anywhere once a
       single token is wider than the column, which is what pre-wrap with
       overflow-wrap:anywhere does. */
    function rowsFor(line, cols) {
        var len = line.length;
        if (len === 0) return 1;
        if (line.indexOf('\t') >= 0) {
            var col = 0;
            for (var t = 0; t < len; t++) {
                col = (line.charAt(t) === '\t') ? (col + 8 - (col % 8)) : col + 1;
            }
            return Math.max(1, Math.ceil(col / cols));
        }
        if (len <= cols) return 1;
        var rows = 1, filled = 0, i = 0;
        while (i < len) {
            var j = line.indexOf(' ', i);
            j = (j < 0) ? len : j + 1;
            var word = j - i;
            if (word > cols) {
                if (filled > 0) { rows++; filled = 0; }
                rows += Math.ceil(word / cols) - 1;
                filled = word % cols;
                if (filled === 0) filled = cols;
            } else if (filled + word > cols) {
                rows++;
                filled = word;
            } else {
                filled += word;
            }
            i = j;
        }
        return rows;
    }

    function blockText(st, k) {
        var b = st.blocks[k];
        var s = st.lines.slice(b.first, b.first + b.count).join('\n');
        /* Every block but the last carries the newline that separated it from
           the next one, so the blocks read back as the original string. */
        return (k < st.blocks.length - 1) ? (s + '\n') : s;
    }

    function applyHeights(el) {
        var st = el.__dtv;
        if (!st) return;
        var m = metrics(el);
        var kids = el.children;
        for (var k = 0; k < st.blocks.length; k++) {
            var b = st.blocks[k];
            var rows = 0;
            for (var i = b.first; i < b.first + b.count; i++) {
                rows += rowsFor(st.lines[i], m.cols);
            }
            b.rows = rows;
            if (kids[k]) {
                kids[k].style.containIntrinsicSize =
                    '0 ' + (rows * m.lineHeight).toFixed(2) + 'px';
            }
        }
    }

    /* The predicted heights are for one column width. Dragging the splitter
       changes it, so recompute -- debounced, because a drag is a stream of
       widths and only the one it stops at matters. */
    function watchResize(el) {
        if (el.__dtvResize || typeof ResizeObserver === 'undefined') return;
        var timer = null;
        var obs = new ResizeObserver(function () {
            if (!el.__dtv) return;
            var m = el.__dtvMetrics;
            if (m && m.width === el.clientWidth) return;
            if (timer) clearTimeout(timer);
            timer = setTimeout(function () {
                timer = null;
                if (el.__dtv) applyHeights(el);
            }, 150);
        });
        obs.observe(el);
        el.__dtvResize = obs;
    }

    function blockIndexAt(st, offset) {
        var starts = st.charStarts;
        var lo = 0, hi = starts.length - 1;
        while (lo < hi) {
            var mid = (lo + hi + 1) >> 1;
            if (starts[mid] <= offset) lo = mid; else hi = mid - 1;
        }
        return lo;
    }

    function restoreBlock(el, k) {
        var st = el.__dtv;
        var node = el.children[k];
        if (!node) return;
        node.textContent = blockText(st, k);
        node.style.contentVisibility = '';
    }

    function clearMarks(el) {
        var st = el && el.__dtv;
        if (!st) return;
        for (var i = 0; i < st.marked.length; i++) restoreBlock(el, st.marked[i]);
        st.marked = [];
    }

    /* Marks go in by rewriting only the blocks that hold one. The block is
       pinned to `visible` while it is marked, because a skipped block lays out
       none of its children and the caller needs the mark's position to scroll
       to it. */
    function markBlock(el, k, ranges, currentId) {
        var st = el.__dtv;
        var node = el.children[k];
        if (!node) return null;
        var text = blockText(st, k);
        var html = '';
        var last = 0;
        var first = null;
        for (var i = 0; i < ranges.length; i++) {
            var s = Math.max(0, Math.min(ranges[i].start, text.length));
            var e = Math.max(s, Math.min(ranges[i].end, text.length));
            if (e <= s) continue;
            html += escapeHtml(text.slice(last, s));
            var cls = 'dtv-mark' + (ranges[i].current ? ' dtv-current' : '');
            var id = '';
            if (ranges[i].current && currentId) {
                /* Both, because the two viewers look the current hit up
                   differently: the Remote Archive by id, and the Calculations
                   browser by class -- its scripts are rewritten to query
                   within one tab's subtree, since three copies of that tab can
                   be on the page and an id would answer for the first. */
                cls += ' ' + currentId;
                id = ' id="' + currentId + '"';
            }
            html += '<mark class="' + cls + '"' + id + '>' + escapeHtml(text.slice(s, e)) + '</mark>';
            last = e;
        }
        html += escapeHtml(text.slice(last));
        node.style.contentVisibility = 'visible';
        node.innerHTML = html;
        if (st.marked.indexOf(k) < 0) st.marked.push(k);
        first = node.querySelector('.dtv-current') || node.querySelector('.dtv-mark');
        return first;
    }

    var TV = {};

    TV.setText = function (el, text, opts) {
        if (!el) return;
        opts = opts || {};
        var lines = String(text == null ? '' : text).split('\n');
        var m = metrics(el);

        var blocks = [];
        var charStarts = [];
        var first = 0, count = 0, chars = 0, rows = 0, offset = 0;
        for (var i = 0; i < lines.length; i++) {
            var line = lines[i];
            count++;
            chars += line.length + 1;
            rows += rowsFor(line, m.cols);
            if (count >= BLOCK_LINES || chars >= BLOCK_CHARS || i === lines.length - 1) {
                blocks.push({first: first, count: count, rows: rows});
                charStarts.push(offset);
                offset += chars;
                first = i + 1;
                count = 0; chars = 0; rows = 0;
            }
        }
        /* The last block does not carry a trailing newline. */
        if (blocks.length) offset -= 1;

        el.__dtv = {
            lines: lines,
            blocks: blocks,
            charStarts: charStarts,
            length: Math.max(0, offset),
            marked: []
        };

        var frag = document.createDocumentFragment();
        for (var k = 0; k < blocks.length; k++) {
            var node = document.createElement('div');
            node.className = 'dtv-blk';
            node.style.containIntrinsicSize =
                '0 ' + (blocks[k].rows * m.lineHeight).toFixed(2) + 'px';
            node.textContent = blockText(el.__dtv, k);
            frag.appendChild(node);
        }
        el.textContent = '';
        el.appendChild(frag);
        watchResize(el);
    };

    /* The text as it went in. Callers used to read el.textContent, which still
       works, but this does not walk the DOM to answer. */
    TV.getText = function (el) {
        var st = el && el.__dtv;
        if (!st) return (el && el.textContent) || '';
        return st.lines.join('\n');
    };

    TV.length = function (el) {
        var st = el && el.__dtv;
        return st ? st.length : ((el && el.textContent) || '').length;
    };

    TV.clearMarks = clearMarks;

    /* Mark one span, given as character offsets into the text that was set.
       Returns the mark element, or null. A span may straddle a block
       boundary, in which case each block gets its own piece of it. */
    TV.markSpan = function (el, start, end, opts) {
        var st = el && el.__dtv;
        if (!st) return null;
        opts = opts || {};
        clearMarks(el);
        var s = Math.max(0, Math.min(start | 0, st.length));
        var e = Math.max(s, Math.min(end | 0, st.length));
        if (e <= s) return null;
        var k0 = blockIndexAt(st, s);
        var k1 = blockIndexAt(st, e - 1);
        var firstMark = null;
        for (var k = k0; k <= k1; k++) {
            var base = st.charStarts[k];
            var mark = markBlock(el, k, [{
                start: Math.max(0, s - base),
                end: e - base,
                current: true
            }], opts.currentId);
            if (!firstMark) firstMark = mark;
        }
        return firstMark;
    };

    /* Every occurrence of a query, for the files small enough that the viewer
       highlights all of them at once. */
    TV.markAll = function (el, query, currentIndex, opts) {
        var st = el && el.__dtv;
        if (!st) return null;
        opts = opts || {};
        clearMarks(el);
        if (!query) return null;
        var text = st.lines.join('\n');
        var re = new RegExp(query.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'), 'gi');
        var byBlock = {};
        var order = [];
        var m, i = 0;
        while ((m = re.exec(text)) !== null) {
            if (m[0].length === 0) { re.lastIndex++; continue; }
            var k = blockIndexAt(st, m.index);
            var base = st.charStarts[k];
            if (!byBlock[k]) { byBlock[k] = []; order.push(k); }
            byBlock[k].push({
                start: m.index - base,
                end: m.index + m[0].length - base,
                current: (i === currentIndex)
            });
            i++;
        }
        var currentMark = null;
        for (var n = 0; n < order.length; n++) {
            var key = order[n];
            var mark = markBlock(el, key, byBlock[key], opts.currentId);
            if (!currentMark && opts.currentId) {
                currentMark = document.getElementById(opts.currentId) || null;
            }
            if (!currentMark) currentMark = mark;
        }
        return opts.currentId
            ? (document.getElementById(opts.currentId) || currentMark)
            : currentMark;
    };

    /* Put a mark in the middle of its scroller. Kept here so both viewers
       centre a hit the same way. */
    TV.centreOn = function (box, mark) {
        if (!box || !mark) return;
        var boxRect = box.getBoundingClientRect();
        var markRect = mark.getBoundingClientRect();
        box.scrollTop += (markRect.top - boxRect.top) - (box.clientHeight / 2);
    };

    window.__delfinTextView = TV;
})();
"""
