#!/usr/bin/env python3
"""Build the README header animation for a DELFIN-family logo.

Lives next to the logos it operates on. Every brand in this folder was
produced with this script -- the only things that change are the logo file,
the wordmark, its colour, and the canvas width needed to fit the letters.

The commands that reproduce the checked-in animations, run from the repo
root, are:

  DELFIN
    python delfin/logo/build_readme_logo_gif.py

  MANTA
    python delfin/logo/build_readme_logo_gif.py \
        --logo delfin/logo/MANTA_LOGO.png --text MANTA \
        --color '#49bbc5' --out delfin/logo/MANTA_readme_demo.gif

  WEDDELL
    python delfin/logo/build_readme_logo_gif.py \
        --logo delfin/logo/WEDDELL_LOGO.png --text WEDDELL \
        --color '#745944' --width 1120 --font-size 108 \
        --out delfin/logo/WEDDELL_readme_demo.gif

  ChemDarwin
    python delfin/logo/build_readme_logo_gif.py \
        --logo delfin/logo/ChemDarwin_Logo.png --text ChemDarwin \
        --color '#309556' --font-size 92 --width 1180 \
        --out delfin/logo/ChemDarwin_readme_demo.gif

For a new brand: pick the wordmark colour from the logo itself so the two
read as one mark, and widen the canvas until the letters are not cramped --
those are the only two values that ever need tuning.

The logo PNG should be cropped flush to its artwork and transparent
outside it, with white underneath the transparent pixels: the frames are
composed on white and downscaled with LANCZOS, which blends whatever RGB
sits under alpha 0 into the edge. A tinted or black background there shows
up as a halo or a dark ring around the logo.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.strip().lstrip("#")
    if len(value) != 6:
        raise ValueError(f"Expected a 6-digit hex color, got {value!r}")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def mix(a: tuple[int, int, int], b: tuple[int, int, int], p: float) -> tuple[int, int, int]:
    p = max(0.0, min(1.0, p))
    return tuple(int(a[i] + (b[i] - a[i]) * p) for i in range(3))


def lighten(color: tuple[int, int, int], amount: float = 0.30) -> tuple[int, int, int]:
    return mix(color, (255, 255, 255), amount)


def darken(color: tuple[int, int, int], amount: float = 0.28) -> tuple[int, int, int]:
    return mix(color, (0, 0, 0), amount)


def ease_out(p: float) -> float:
    p = max(0.0, min(1.0, p))
    return 1 - (1 - p) ** 3


def ease_in_out(p: float) -> float:
    p = max(0.0, min(1.0, p))
    return 0.5 - 0.5 * math.cos(math.pi * p)


def load_font(size: int) -> ImageFont.FreeTypeFont:
    """Return DejaVu Sans Bold, the face every checked-in animation uses.

    The distributions disagree on where it lives, so all known layouts are
    tried, then matplotlib's bundled copy. Falling back to Pillow's default
    is not acceptable here: it is a 10 px bitmap face, so the wordmark would
    come out unreadably small instead of failing visibly. Better to stop.
    """
    candidates = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",   # Debian, Ubuntu
        "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",            # Arch
        "/usr/share/fonts/dejavu-sans-fonts/DejaVuSans-Bold.ttf",  # RHEL, Rocky
        "/usr/share/fonts/TTF/DejaVuSans-Bold.ttf",
        "/Library/Fonts/DejaVuSans-Bold.ttf",                     # macOS
    ]
    try:
        import matplotlib
        candidates.append(
            str(Path(matplotlib.get_data_path()) / "fonts" / "ttf" / "DejaVuSans-Bold.ttf")
        )
    except Exception:
        pass

    for path in candidates:
        try:
            return ImageFont.truetype(path, size)
        except OSError:
            continue
    # Pillow also searches the system font directories by bare filename.
    try:
        return ImageFont.truetype("DejaVuSans-Bold.ttf", size)
    except OSError:
        raise SystemExit(
            "DejaVuSans-Bold.ttf not found. Install the DejaVu fonts "
            "(dejavu-sans-fonts / fonts-dejavu) or pip install matplotlib, "
            "which ships a copy."
        )


def build_gif(
    *,
    logo_path: Path,
    text: str,
    out_path: Path,
    color: tuple[int, int, int],
    width: int,
    height: int,
    font_size: int,
    left_shift: int,
) -> None:
    scale_factor = 4
    sw, sh = width * scale_factor, height * scale_factor
    s_left_shift = left_shift * scale_factor

    logo = Image.open(logo_path).convert("RGBA")
    font = load_font(font_size * scale_factor)
    light_color = lighten(color, 0.32)
    shadow_color = darken(color, 0.34)

    spacing = 8 * scale_factor if len(text) <= 6 else 4 * scale_factor
    probe = ImageDraw.Draw(Image.new("RGB", (1, 1), "white"))
    widths: list[int] = []
    for ch in text:
        box = probe.textbbox((0, 0), ch, font=font)
        widths.append(box[2] - box[0])

    text_width = sum(widths) + spacing * (len(text) - 1)
    logo_size = 172 * scale_factor
    gap = 42 * scale_factor
    group_width = logo_size + gap + text_width
    group_left = (sw - group_width) // 2 - s_left_shift
    logo_center = (group_left + logo_size // 2, sh // 2)
    text_x = group_left + logo_size + gap
    text_y = (68 if font_size >= 110 else 82) * scale_factor

    frames_rgb: list[Image.Image] = []
    frame_count = 58 if len(text) <= 6 else 66
    text_step = 0.068 if len(text) <= 6 else 0.045
    settle = 0.17 if len(text) <= 6 else 0.15

    for frame in range(frame_count):
        t = frame / (frame_count - 1)
        image = Image.new("RGB", (sw, sh), "white")

        logo_progress = ease_out(t / 0.20)
        logo_scale = 0.96 + 0.04 * logo_progress + 0.003 * math.sin(t * math.tau)
        logo_frame = logo.copy()
        target = int(logo_size * logo_scale)
        logo_frame = logo_frame.resize(
            (target, int(target * logo_frame.height / logo_frame.width)),
            Image.Resampling.LANCZOS,
        )

        rgba = image.convert("RGBA")
        rgba.paste(
            logo_frame,
            (logo_center[0] - logo_frame.width // 2, logo_center[1] - logo_frame.height // 2),
            logo_frame,
        )
        image = rgba.convert("RGB")

        x = text_x
        for idx, ch in enumerate(text):
            start = 0.20 + idx * text_step
            p = ease_out((t - start) / settle)
            if p <= 0:
                x += widths[idx] + spacing
                continue

            lift = math.sin(min(1.0, p) * math.pi)
            y = text_y + int((1 - p) * 26 * scale_factor - lift * 8 * scale_factor)
            letter_color = mix(light_color, color, ease_in_out((p - 0.28) / 0.72))

            shadow = Image.new("RGBA", (sw, sh), (0, 0, 0, 0))
            shadow_draw = ImageDraw.Draw(shadow)
            shadow_alpha = int(36 * p * (1.0 - 0.35 * lift))
            shadow_draw.ellipse(
                (
                    x + 8 * scale_factor,
                    text_y + 118 * scale_factor,
                    x + widths[idx] - 7 * scale_factor,
                    text_y + 133 * scale_factor,
                ),
                fill=shadow_color + (shadow_alpha,),
            )
            shadow = shadow.filter(ImageFilter.GaussianBlur(7 * scale_factor))
            image = Image.alpha_composite(image.convert("RGBA"), shadow).convert("RGB")

            ImageDraw.Draw(image).text((x, y), ch, font=font, fill=letter_color)
            x += widths[idx] + spacing

        frames_rgb.append(image.resize((width, height), Image.Resampling.LANCZOS))

    frames_rgb.extend([frames_rgb[-1]] * 16)

    rep_indices = sorted(
        set([0, 4, 8, 12, 16, 22, 30, 40, frame_count - 1, len(frames_rgb) - 1])
    )
    palette_source = Image.new("RGB", (width, height * len(rep_indices)), "white")
    for row, idx in enumerate(rep_indices):
        palette_source.paste(frames_rgb[idx], (0, height * row))
    palette = palette_source.quantize(
        colors=256,
        method=Image.Quantize.MEDIANCUT,
        dither=Image.Dither.NONE,
    )
    frames = [frame.quantize(palette=palette, dither=Image.Dither.NONE) for frame in frames_rgb]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(
        out_path,
        save_all=True,
        append_images=frames[1:],
        duration=64,
        loop=0,
        optimize=True,
    )


def main() -> None:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--logo", type=Path, default=here / "DELFIN_logo.png")
    parser.add_argument("--text", default="DELFIN")
    parser.add_argument("--out", type=Path, default=here / "DELFIN_readme_demo.gif")
    parser.add_argument("--color", default="#628ebe",
                        help="wordmark colour; match it to the logo")
    parser.add_argument("--width", type=int, default=1050)
    parser.add_argument("--height", type=int, default=280)
    parser.add_argument("--font-size", type=int, default=116)
    parser.add_argument("--left-shift", type=int, default=15)
    args = parser.parse_args()

    if not args.logo.is_file():
        raise SystemExit(f"logo not found: {args.logo}")

    build_gif(
        logo_path=args.logo,
        text=args.text,
        out_path=args.out,
        color=hex_to_rgb(args.color),
        width=args.width,
        height=args.height,
        font_size=args.font_size,
        left_shift=args.left_shift,
    )
    print(args.out)


if __name__ == "__main__":
    sys.exit(main())
