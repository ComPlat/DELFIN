"""A QR code encoder for one purpose: the phone-pairing link, as a picture.

Byte mode, error correction level L or M, versions 1-40, the eight masks
scored by the standard's penalty rules (ISO/IEC 18004). Pure Python, no
dependency, so every host can draw the code (qrencode where installed is
not needed). Input: a text; output: a square matrix of booleans (True is a
dark module), without the quiet zone.
"""
from __future__ import annotations

# Error correction codewords per block and number of blocks, per version
# (index 0 = version 1), levels L and M (ISO/IEC 18004 table 9).
_ECC_PER_BLOCK = {
    "L": [7, 10, 15, 20, 26, 18, 20, 24, 30, 18, 20, 24, 26, 30, 22, 24, 28, 30, 28, 28,
          28, 28, 30, 30, 26, 28, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30],
    "M": [10, 16, 26, 18, 24, 16, 18, 22, 22, 26, 30, 22, 22, 24, 24, 28, 28, 26, 26, 26,
          26, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28],
}
_NUM_BLOCKS = {
    "L": [1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 7, 8,
          8, 9, 9, 10, 12, 12, 12, 13, 14, 15, 16, 17, 18, 19, 19, 20, 21, 22, 24, 25],
    "M": [1, 1, 1, 2, 2, 4, 4, 4, 5, 5, 5, 8, 9, 9, 10, 10, 11, 13, 14, 16,
          17, 17, 18, 20, 21, 23, 25, 26, 28, 29, 31, 33, 35, 37, 38, 40, 43, 45, 47, 49],
}


def _raw_modules(version: int) -> int:
    """Modules left for data and error correction in a version."""
    result = (16 * version + 128) * version + 64
    if version >= 2:
        n = version // 7 + 2
        result -= (25 * n - 10) * n - 55
        if version >= 7:
            result -= 36
    return result


_ALIGN = [
    [], [6, 18], [6, 22], [6, 26], [6, 30], [6, 34], [6, 22, 38], [6, 24, 42], [6, 26, 46], [6, 28, 50],
    [6, 30, 54], [6, 32, 58], [6, 34, 62], [6, 26, 46, 66], [6, 26, 48, 70], [6, 26, 50, 74], [6, 30, 54, 78],
    [6, 30, 56, 82], [6, 30, 58, 86], [6, 34, 62, 90], [6, 28, 50, 72, 94], [6, 26, 50, 74, 98],
    [6, 30, 54, 78, 102], [6, 28, 54, 80, 106], [6, 32, 58, 84, 110], [6, 30, 58, 86, 114], [6, 34, 62, 90, 118],
    [6, 26, 50, 74, 98, 122], [6, 30, 54, 78, 102, 126], [6, 26, 52, 78, 104, 130], [6, 30, 56, 82, 108, 134],
    [6, 34, 60, 86, 112, 138], [6, 30, 58, 86, 114, 142], [6, 34, 62, 90, 118, 146], [6, 30, 54, 78, 102, 126, 150],
    [6, 24, 50, 76, 102, 128, 154], [6, 28, 54, 80, 106, 132, 158], [6, 32, 58, 84, 110, 136, 162],
    [6, 26, 54, 82, 110, 138, 166], [6, 30, 58, 86, 114, 142, 170],
]

# GF(256) with the QR polynomial 0x11d.
_EXP = [0] * 512
_LOG = [0] * 256
_x = 1
for _i in range(255):
    _EXP[_i] = _x
    _LOG[_x] = _i
    _x <<= 1
    if _x & 0x100:
        _x ^= 0x11D
for _i in range(255, 512):
    _EXP[_i] = _EXP[_i - 255]


def _gf_mul(a: int, b: int) -> int:
    if a == 0 or b == 0:
        return 0
    return _EXP[_LOG[a] + _LOG[b]]


def _generator(n: int) -> list[int]:
    poly = [1]
    for i in range(n):
        nxt = [0] * (len(poly) + 1)
        for j, c in enumerate(poly):
            nxt[j] ^= c
            nxt[j + 1] ^= _gf_mul(c, _EXP[i])
        poly = nxt
    return poly


def _ec_codewords(data: list[int], n: int) -> list[int]:
    gen = _generator(n)
    rem = list(data) + [0] * n
    for i in range(len(data)):
        coef = rem[i]
        if coef:
            for j in range(1, len(gen)):
                rem[i + j] ^= _gf_mul(gen[j], coef)
    return rem[len(data):]


def _blocks(level: str, version: int):
    """(data codewords, ec codewords per block, data codewords of each block)."""
    raw = _raw_modules(version) // 8
    ec = _ECC_PER_BLOCK[level][version - 1]
    n = _NUM_BLOCKS[level][version - 1]
    short = n - raw % n
    short_len = raw // n
    sizes = [short_len - ec] * short + [short_len + 1 - ec] * (n - short)
    return sum(sizes), ec, sizes


def _char_count_bits(version: int) -> int:
    return 8 if version <= 9 else 16


def _choose_version(nbytes: int, level: str) -> int:
    for v in range(1, 41):
        data_total, _, _ = _blocks(level, v)
        if 4 + _char_count_bits(v) + 8 * nbytes <= 8 * data_total:
            return v
    raise ValueError("text too long for a QR code")


def _codewords(payload: bytes, level: str, version: int) -> list[int]:
    data_total, ec, sizes = _blocks(level, version)
    bits = [0, 1, 0, 0]                                   # byte mode
    n = len(payload)
    bits += [(n >> i) & 1 for i in range(_char_count_bits(version) - 1, -1, -1)]
    for b in payload:
        bits += [(b >> i) & 1 for i in range(7, -1, -1)]
    cap = 8 * data_total
    bits += [0] * min(4, cap - len(bits))
    bits += [0] * (-len(bits) % 8)
    data = [int("".join(map(str, bits[i:i + 8])), 2) for i in range(0, len(bits), 8)]
    pad = [0xEC, 0x11]
    i = 0
    while len(data) < data_total:
        data.append(pad[i % 2])
        i += 1
    blocks, pos = [], 0
    for size in sizes:
        blocks.append(data[pos:pos + size])
        pos += size
    ecs = [_ec_codewords(b, ec) for b in blocks]
    out = []
    for i in range(max(sizes)):
        for b in blocks:
            if i < len(b):
                out.append(b[i])
    for i in range(ec):
        for e in ecs:
            out.append(e[i])
    return out


def _remainder_bits(version: int) -> int:
    if version in (2, 3, 4, 5, 6):
        return 7
    if 14 <= version <= 20 or 28 <= version <= 34:
        return 3
    if 21 <= version <= 27:
        return 4
    return 0


def _bch_format(level: str, mask: int) -> int:
    data = ({"L": 1, "M": 0}[level] << 3) | mask
    rem = data << 10
    for i in range(14, 9, -1):
        if rem & (1 << i):
            rem ^= 0x537 << (i - 10)
    return ((data << 10) | rem) ^ 0x5412


def _bch_version(version: int) -> int:
    rem = version << 12
    for i in range(17, 11, -1):
        if rem & (1 << i):
            rem ^= 0x1F25 << (i - 12)
    return (version << 12) | rem


_MASKS = [
    lambda r, c: (r + c) % 2 == 0,
    lambda r, c: r % 2 == 0,
    lambda r, c: c % 3 == 0,
    lambda r, c: (r + c) % 3 == 0,
    lambda r, c: (r // 2 + c // 3) % 2 == 0,
    lambda r, c: (r * c) % 2 + (r * c) % 3 == 0,
    lambda r, c: ((r * c) % 2 + (r * c) % 3) % 2 == 0,
    lambda r, c: ((r + c) % 2 + (r * c) % 3) % 2 == 0,
]


def _base(version: int):
    size = 17 + 4 * version
    m = [[None] * size for _ in range(size)]

    def finder(r0, c0):
        for r in range(-1, 8):
            for c in range(-1, 8):
                rr, cc = r0 + r, c0 + c
                if 0 <= rr < size and 0 <= cc < size:
                    dark = (0 <= r <= 6 and c in (0, 6)) or (0 <= c <= 6 and r in (0, 6)) \
                        or (2 <= r <= 4 and 2 <= c <= 4)
                    m[rr][cc] = dark
    finder(0, 0)
    finder(0, size - 7)
    finder(size - 7, 0)
    for i in range(8, size - 8):
        m[6][i] = m[i][6] = i % 2 == 0
    pos = _ALIGN[version - 1]
    last = len(pos) - 1
    for i, r in enumerate(pos):
        for j, c in enumerate(pos):
            if (i, j) in {(0, 0), (0, last), (last, 0)}:
                continue                                  # the finder corners
            for dr in range(-2, 3):
                for dc in range(-2, 3):
                    m[r + dr][c + dc] = max(abs(dr), abs(dc)) != 1
    m[size - 8][8] = True                                  # the dark module
    # format and version areas reserved
    for i in range(9):
        if m[8][i] is None:
            m[8][i] = False
        if m[i][8] is None:
            m[i][8] = False
    for i in range(8):
        if m[8][size - 1 - i] is None:
            m[8][size - 1 - i] = False
        if m[size - 1 - i][8] is None:
            m[size - 1 - i][8] = False
    if version >= 7:
        for i in range(6):
            for j in range(3):
                m[i][size - 11 + j] = False
                m[size - 11 + j][i] = False
    return m


def _place(m, bits):
    size = len(m)
    reserved = [[v is not None for v in row] for row in m]
    i = 0
    col = size - 1
    up = True
    while col > 0:
        if col == 6:
            col -= 1
        rows = range(size - 1, -1, -1) if up else range(size)
        for r in rows:
            for c in (col, col - 1):
                if not reserved[r][c]:
                    m[r][c] = bits[i] if i < len(bits) else False
                    i += 1
        up = not up
        col -= 2
    return reserved


def _apply_format(m, level, mask, version):
    size = len(m)
    f = _bch_format(level, mask)
    bit = [(f >> i) & 1 == 1 for i in range(15)]
    for i in range(6):
        m[i][8] = bit[i]
    m[7][8] = bit[6]
    m[8][8] = bit[7]
    m[8][7] = bit[8]
    for i in range(9, 15):
        m[8][14 - i] = bit[i]
    for i in range(8):
        m[8][size - 1 - i] = bit[i]
    for i in range(8, 15):
        m[size - 15 + i][8] = bit[i]
    m[size - 8][8] = True
    if version >= 7:
        v = _bch_version(version)
        for i in range(18):
            bit = (v >> i) & 1 == 1
            a, b = i // 3, size - 11 + i % 3
            m[a][b] = bit
            m[b][a] = bit


def _penalty(m) -> int:
    size = len(m)
    score = 0
    for grid in (m, [list(col) for col in zip(*m)]):
        for row in grid:
            run, last = 0, None
            for v in row:
                if v == last:
                    run += 1
                else:
                    if run >= 5:
                        score += run - 2
                    run, last = 1, v
            if run >= 5:
                score += run - 2
            line = "".join("1" if v else "0" for v in row)
            for pat in ("10111010000", "00001011101"):
                score += 40 * sum(1 for i in range(len(line) - 10) if line.startswith(pat, i))
    for r in range(size - 1):
        for c in range(size - 1):
            v = m[r][c]
            if v == m[r][c + 1] == m[r + 1][c] == m[r + 1][c + 1]:
                score += 3
    dark = sum(v for row in m for v in row)
    score += 10 * (abs(dark * 20 - size * size * 10) // (size * size))
    return score


def encode(text: str, level: str = "L") -> list[list[bool]]:
    """``text`` (UTF-8, byte mode) as a QR matrix, the smallest version that
    holds it, the mask with the lowest penalty."""
    payload = text.encode("utf-8")
    version = _choose_version(len(payload), level)
    words = _codewords(payload, level, version)
    bits = [(w >> i) & 1 == 1 for w in words for i in range(7, -1, -1)]
    bits += [False] * _remainder_bits(version)
    best = None
    for mask in range(8):
        m = _base(version)
        reserved = _place(m, bits)
        rule = _MASKS[mask]
        for r in range(len(m)):
            for c in range(len(m)):
                if not reserved[r][c] and rule(r, c):
                    m[r][c] = not m[r][c]
        _apply_format(m, level, mask, version)
        score = _penalty(m)
        if best is None or score < best[0]:
            best = (score, m)
    return [[bool(v) for v in row] for row in best[1]]


def svg(text: str, *, scale: int = 6, border: int = 4) -> str:
    """The code as an SVG image (dark on white, a quiet zone of ``border``)."""
    m = encode(text)
    n = len(m) + 2 * border
    rects = "".join(f'M{c + border},{r + border}h1v1h-1z'
                    for r, row in enumerate(m) for c, v in enumerate(row) if v)
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {n} {n}" '
            f'width="{n * scale}" height="{n * scale}" shape-rendering="crispEdges">'
            f'<rect width="{n}" height="{n}" fill="#fff"/><path d="{rects}" fill="#000"/></svg>')


def text_blocks(text: str, *, border: int = 2) -> str:
    """The code as text with half blocks (two rows per line), for a terminal."""
    m = encode(text)
    n = len(m)
    grid = [[False] * (n + 2 * border) for _ in range(n + 2 * border)]
    for r in range(n):
        for c in range(n):
            grid[r + border][c + border] = m[r][c]
    if len(grid) % 2:
        grid.append([False] * len(grid[0]))
    lines = []
    for r in range(0, len(grid), 2):
        line = []
        for c in range(len(grid[0])):
            top, bottom = grid[r][c], grid[r + 1][c]
            # Light text on a dark terminal shows dark modules as spaces:
            # draw the LIGHT modules, so the code reads right on dark themes.
            line.append({(False, False): "█", (True, False): "▄", (False, True): "▀", (True, True): " "}[(top, bottom)])
        lines.append("".join(line))
    return "\n".join(lines)
