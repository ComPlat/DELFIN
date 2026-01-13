"""ESP report helpers.

Primary target: a publication-style 3D rendering showing the molecule plus an electron-density
isosurface colored by the electrostatic potential (ESP), with a colorbar embedded into a single PNG.

Fallback: If PyMOL or an electron-density cube is unavailable, we fall back to a 2D mid-plane slice
of the ESP cube with matplotlib.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Optional, Sequence

from delfin.common.logging import get_logger

logger = get_logger(__name__)

BOHR_TO_ANGSTROM = 0.529177210903
DEFAULT_ESP_VLIM_ABS = 0.025  # a.u.; strong contrast in reports


@dataclass(frozen=True)
class CubeHeader:
    natoms: int
    origin_bohr: tuple[float, float, float]
    nx: int
    ny: int
    nz: int
    vx_bohr: tuple[float, float, float]
    vy_bohr: tuple[float, float, float]
    vz_bohr: tuple[float, float, float]


def _find_first_existing(workspace: Path, candidates: list[str]) -> Optional[Path]:
    for name in candidates:
        for folder in (workspace / "ESD", workspace):
            candidate = folder / name
            if candidate.exists():
                return candidate
    return None


def _orca_plot_binary() -> Optional[str]:
    """Return a usable orca_plot executable if present."""
    direct = shutil.which("orca_plot")
    if direct:
        return direct
    fallback = Path("/opt/orca/orca_plot")
    if fallback.exists():
        return str(fallback)
    return None


def _orca_plot_run(gbw_path: Path, stdin: str) -> bool:
    """Run orca_plot in interactive/batch mode with a provided stdin script."""
    exe = _orca_plot_binary()
    if not exe:
        logger.warning("orca_plot not found (PATH or /opt/orca/orca_plot); skipping cube generation")
        return False

    if not gbw_path.exists():
        return False

    cmd = [exe, str(gbw_path.name), "-i"]
    logger.info("Running %s (cwd=%s)", " ".join(cmd), gbw_path.parent)
    try:
        start = time.time()
        proc = subprocess.run(
            cmd,
            input=stdin,
            text=True,
            capture_output=True,
            cwd=str(gbw_path.parent),
            check=False,
        )
        elapsed = time.time() - start
        if proc.returncode != 0:
            logger.error(
                "orca_plot failed (rc=%s, %.1fs).\nstdout:\n%s\nstderr:\n%s",
                proc.returncode,
                elapsed,
                (proc.stdout or "").strip(),
                (proc.stderr or "").strip(),
            )
            return False
        logger.info("orca_plot finished in %.1fs", elapsed)
        return True
    except Exception as exc:  # noqa: BLE001
        logger.error("orca_plot execution failed: %s", exc, exc_info=True)
        return False


def _orca_plot_generate_esp_cube(gbw_path: Path) -> bool:
    """Generate an ESP cube via orca_plot.

    Sequence (as requested):
      1
      43
      0
      11
      12
    """
    return _orca_plot_run(gbw_path, "1\n43\n0\n11\n12\n")


def _pick_esp_cube(folder: Path, modified_after: float | None = None) -> Optional[Path]:
    """Pick the most relevant ESP cube in a folder."""
    candidates = list(folder.glob("*.esp.cube"))
    if not candidates:
        candidates = list(folder.glob("*0.esp.cube")) + list(folder.glob("0.esp.cube"))
    if modified_after is not None:
        candidates = [p for p in candidates if p.stat().st_mtime >= modified_after]
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def _pick_density_cube(folder: Path, modified_after: float | None = None) -> Optional[Path]:
    """Pick a likely electron-density cube in a folder."""
    patterns = [
        "*eldens*.cube",
        "*edens*.cube",
        "*density*.cube",
        "*dens*.cube",
        "*rho*.cube",
    ]
    candidates: list[Path] = []
    for pat in patterns:
        candidates.extend(folder.glob(pat))

    # Filter out common non-density cubes
    filtered: list[Path] = []
    for p in candidates:
        name = p.name.lower()
        if "esp" in name:
            continue
        if ".mo" in name or "mo" in name:
            continue
        if "spin" in name or "sden" in name:
            continue
        filtered.append(p)

    if modified_after is not None:
        filtered = [p for p in filtered if p.stat().st_mtime >= modified_after]
    if not filtered:
        return None
    return max(filtered, key=lambda p: p.stat().st_mtime)


def _orca_plot_generate_density_cube(gbw_path: Path, plot_type_ids: Sequence[int] = (41, 40, 42, 44)) -> bool:
    """Try to generate an electron-density cube.

    ORCA's plot-type numbering can differ between versions; we try a small set of common IDs.
    """
    for plot_type in plot_type_ids:
        ok = _orca_plot_run(gbw_path, f"1\n{plot_type}\n0\n11\n12\n")
        if ok:
            return True
    return False


def _parse_cube_header(lines: list[str]) -> CubeHeader:
    # Cube format: 2 comment lines, then:
    # natoms ox oy oz
    # nx vx1 vx2 vx3
    # ny vy1 vy2 vy3
    # nz vz1 vz2 vz3
    third = lines[2].split()
    natoms = int(float(third[0]))
    origin = tuple(float(x) for x in third[1:4])

    nx_line = lines[3].split()
    ny_line = lines[4].split()
    nz_line = lines[5].split()

    nx = int(float(nx_line[0]))
    ny = int(float(ny_line[0]))
    nz = int(float(nz_line[0]))

    vx = tuple(float(x) for x in nx_line[1:4])
    vy = tuple(float(x) for x in ny_line[1:4])
    vz = tuple(float(x) for x in nz_line[1:4])

    return CubeHeader(
        natoms=abs(natoms),
        origin_bohr=origin,
        nx=abs(nx),
        ny=abs(ny),
        nz=abs(nz),
        vx_bohr=vx,
        vy_bohr=vy,
        vz_bohr=vz,
    )


def _read_cube(cube_path: Path):
    """Read Gaussian cube file as (header, data[z,y,x])."""
    import numpy as np

    text = cube_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if len(text) < 10:
        raise ValueError(f"Cube file too short: {cube_path}")

    header = _parse_cube_header(text[:6])
    data_start = 6 + header.natoms
    if data_start >= len(text):
        raise ValueError(f"Cube file missing data section: {cube_path}")

    raw = " ".join(text[data_start:])
    values = np.fromstring(raw, sep=" ")
    expected = header.nx * header.ny * header.nz
    if values.size < expected:
        raise ValueError(f"Cube data truncated: expected {expected} values, got {values.size}")
    if values.size > expected:
        values = values[:expected]

    # Cube ordering: x fastest, then y, then z.
    data = values.reshape((header.nz, header.ny, header.nx))
    return header, data


def _robust_symmetric_vlim(values) -> float:
    import numpy as np

    finite = np.isfinite(values)
    if not finite.any():
        return 1.0
    abs_v = np.abs(values[finite])
    if abs_v.size == 0:
        return 1.0
    vlim = float(np.nanpercentile(abs_v, 99.0))
    if not vlim or vlim <= 0:
        vlim = float(np.nanmax(abs_v))
    return vlim if vlim and vlim > 0 else 1.0


def _esp_vlim_from_cube(
    esp_values,
    percentile: float = 97.0,
    factor: float = 0.6,
    min_vlim: float = 1e-4,
) -> float:
    """Choose a symmetric color range for ESP that yields visually meaningful contrast.

    Using the raw max often produces a nearly-white surface; a robust percentile with a factor
    gives more saturated red/blue regions.
    """
    import numpy as np

    finite = np.isfinite(esp_values)
    if not finite.any():
        return 0.05
    abs_v = np.abs(esp_values[finite])
    if abs_v.size == 0:
        return 0.05

    try:
        base = float(np.nanpercentile(abs_v, percentile))
    except Exception:
        base = float(np.nanmax(abs_v))

    if not base or base <= 0:
        base = float(np.nanmax(abs_v)) if abs_v.size else 0.05

    vlim = base * float(factor)
    if vlim < min_vlim:
        vlim = max(min_vlim, float(base) if base else min_vlim)
    return float(vlim)


def _create_colorbar_png(output_png: Path, vlim: float, label: str = "Electrostatic potential (a.u.)") -> Optional[Path]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib as mpl
    except Exception as exc:  # noqa: BLE001
        logger.warning("matplotlib not available; cannot render ESP colorbar: %s", exc)
        return None

    # Keep this compact: only the color strip + tick labels (no extra frame).
    fig = plt.figure(figsize=(0.9, 4.4))
    ax = fig.add_axes([0.55, 0.04, 0.22, 0.92])  # [left,bottom,width,height]
    norm = mpl.colors.Normalize(vmin=-vlim, vmax=vlim)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap="coolwarm")
    cb = fig.colorbar(sm, cax=ax, orientation="vertical")
    cb.ax.tick_params(labelsize=8, length=0)
    try:
        cb.outline.set_visible(False)
    except Exception:
        pass
    for spine in cb.ax.spines.values():
        spine.set_visible(False)
    if label:
        # Keep label compact to avoid extra whitespace in tight layouts.
        cb.ax.set_title(label, fontsize=8, pad=2)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_png), dpi=300, bbox_inches="tight", pad_inches=0.0, facecolor="white")
    plt.close(fig)
    return output_png if output_png.exists() else None


def _autocrop_white_pil(im, threshold: int = 10, pad: int = 6):
    try:
        from PIL import Image, ImageChops
    except Exception as exc:  # noqa: BLE001
        logger.warning("Pillow not available; cannot auto-crop whitespace: %s", exc)
        return im

    bg = Image.new("RGB", im.size, "white")
    diff = ImageChops.difference(im, bg).convert("L")
    diff = diff.point(lambda p: 255 if p > threshold else 0)
    bbox = diff.getbbox()
    if not bbox:
        return im
    left, upper, right, lower = bbox
    left = max(0, left - pad)
    upper = max(0, upper - pad)
    right = min(im.width, right + pad)
    lower = min(im.height, lower + pad)
    return im.crop((left, upper, right, lower))


def _merge_with_colorbar(main_png: Path, colorbar_png: Path, output_png: Path, pad_px: int = 8) -> Optional[Path]:
    try:
        from PIL import Image
    except Exception as exc:  # noqa: BLE001
        logger.warning("Pillow not available; cannot merge ESP render + colorbar: %s", exc)
        return main_png

    img = Image.open(main_png).convert("RGB")
    bar = Image.open(colorbar_png).convert("RGB")

    # Remove large white borders (common with ray renders); we will re-add only the right-side colorbar.
    img = _autocrop_white_pil(img, threshold=10, pad=8)
    bar = _autocrop_white_pil(bar, threshold=10, pad=4)

    # Resize bar to match image height
    if bar.height != img.height:
        scale = img.height / float(bar.height)
        new_w = max(1, int(bar.width * scale))
        bar = bar.resize((new_w, img.height), Image.Resampling.LANCZOS)

    pad_px = max(4, int(pad_px))
    out = Image.new("RGB", (img.width + pad_px + bar.width, img.height), "white")
    out.paste(img, (0, 0))
    out.paste(bar, (img.width + pad_px, 0))
    output_png.parent.mkdir(parents=True, exist_ok=True)
    out.save(output_png)
    return output_png if output_png.exists() else None


def _make_grid_png(input_pngs: Sequence[Path], output_png: Path, cols: int = 2, pad_px: int = 10) -> Optional[Path]:
    try:
        from PIL import Image
    except Exception as exc:  # noqa: BLE001
        logger.warning("Pillow not available; cannot create ESP montage: %s", exc)
        return None

    images = []
    for p in input_pngs:
        if not p.exists():
            continue
        img = Image.open(p).convert("RGB")
        img = _autocrop_white_pil(img, threshold=10, pad=8)
        images.append(img)

    if not images:
        return None

    cols = max(1, int(cols))
    rows = (len(images) + cols - 1) // cols
    pad_px = max(0, int(pad_px))

    # Normalize tile sizes to the smallest to avoid upscaling artifacts.
    tile_w = min(im.width for im in images)
    tile_h = min(im.height for im in images)
    images = [im.resize((tile_w, tile_h), Image.Resampling.LANCZOS) for im in images]

    out_w = cols * tile_w + (cols - 1) * pad_px
    out_h = rows * tile_h + (rows - 1) * pad_px
    out = Image.new("RGB", (out_w, out_h), "white")

    for idx, im in enumerate(images):
        r = idx // cols
        c = idx % cols
        x = c * (tile_w + pad_px)
        y = r * (tile_h + pad_px)
        out.paste(im, (x, y))

    output_png.parent.mkdir(parents=True, exist_ok=True)
    out.save(output_png)
    return output_png if output_png.exists() else None


def create_esp_isosurface_png(
    xyz_path: Path,
    density_cube_path: Path,
    esp_cube_path: Path,
    output_png: Path,
    density_isovalue: float = 0.001,
) -> Optional[Path]:
    """Render molecule + density isosurface colored by ESP using PyMOL.

    PyMOL can crash the hosting interpreter on some systems (OpenGL/headless issues).
    To keep DELFIN robust, the actual PyMOL render is executed in a subprocess.
    """

    import os

    try:
        _, esp_data = _read_cube(esp_cube_path)
        vlim_abs_env = os.getenv("DELFIN_ESP_VLIM_ABS")
        if vlim_abs_env:
            vlim = abs(float(vlim_abs_env))
        else:
            vlim = DEFAULT_ESP_VLIM_ABS
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to read ESP cube for scaling (%s): %s", esp_cube_path, exc)
        vlim = DEFAULT_ESP_VLIM_ABS

    output_png.parent.mkdir(parents=True, exist_ok=True)
    # Main render via PyMOL (no GUI) in a subprocess for safety.
    code = r"""
import sys
xyz, dens_cube, esp_cube, out_png, iso_str, vlim_str = sys.argv[1:7]
iso = float(iso_str)
vlim = float(vlim_str)
try:
    import pymol
    from pymol import cmd
except Exception as exc:
    print(f"PyMOL import failed: {exc}", file=sys.stderr)
    sys.exit(2)

try:
    pymol.finish_launching(["pymol", "-cq"])
    cmd.reinitialize()

    cmd.load(xyz, "mol")
    cmd.hide("everything", "mol")
    cmd.show("sticks", "mol")
    cmd.show("spheres", "mol")
    cmd.set("sphere_scale", 0.25)
    cmd.set("stick_radius", 0.15)
    cmd.set("stick_quality", 15)
    cmd.set("sphere_quality", 3)
    cmd.set("valence", 1)

    cmd.color("gray50", "elem C")
    cmd.color("white", "elem H")
    cmd.color("slate", "elem N")
    cmd.color("red", "elem O")
    cmd.color("yellow", "elem S")

    cmd.bg_color("white")
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("ambient", 0.35)
    cmd.set("direct", 0.65)
    cmd.set("specular", 0.2)
    cmd.set("shininess", 10)
    cmd.set("depth_cue", 0)

    cmd.load(dens_cube, "density")
    cmd.load(esp_cube, "esp")

    cmd.isosurface("dens_surf", "density", iso)
    cmd.set("transparency", 0.35, "dens_surf")

    cmd.ramp_new("esp_ramp", "esp", [-vlim, 0.0, vlim], ["red", "white", "blue"])
    cmd.set("surface_color", "esp_ramp", "dens_surf")

    cmd.orient("mol")
    cmd.turn("y", -35)
    cmd.turn("x", 15)
    cmd.zoom("all", buffer=2.5)

    cmd.ray(1800, 1400)
    cmd.png(out_png, dpi=300)
except Exception as exc:
    print(f"PyMOL render failed: {exc}", file=sys.stderr)
    sys.exit(3)
finally:
    try:
        cmd.delete("all")
        cmd.reinitialize()
    except Exception:
        pass
"""
    proc = subprocess.run(
        [sys.executable, "-c", code, str(xyz_path), str(density_cube_path), str(esp_cube_path), str(output_png), str(density_isovalue), str(vlim)],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0 or not output_png.exists():
        logger.warning(
            "PyMOL subprocess render failed (rc=%s).\nstdout:\n%s\nstderr:\n%s",
            proc.returncode,
            (proc.stdout or "").strip(),
            (proc.stderr or "").strip(),
        )
        return None

    if not output_png.exists():
        return None

    import os

    # Add colorbar (optional)
    if str(os.getenv("DELFIN_ESP_ADD_COLORBAR", "1")).lower() in {"0", "false", "no"}:
        return output_png

    colorbar_png = output_png.with_name(output_png.stem + "_colorbar.png")
    cb = _create_colorbar_png(colorbar_png, vlim=vlim)
    if cb and cb.exists():
        merged = _merge_with_colorbar(output_png, cb, output_png)
        try:
            colorbar_png.unlink()
        except Exception:
            pass
        return merged

    return output_png


def create_esp_molecular_surface_png(
    xyz_path: Path,
    esp_cube_path: Path,
    output_png: Path,
    surface_transparency: float = 0.10,
) -> Optional[Path]:
    """Render molecule surface colored by ESP (no density cube required).

    This yields the typical "MEP mapped on molecular surface" style visualization, while staying
    robust even when ORCA does not provide a usable electron-density cube.
    """
    try:
        _, esp_data = _read_cube(esp_cube_path)
        import os

        vlim_abs_env = os.getenv("DELFIN_ESP_VLIM_ABS")
        if vlim_abs_env:
            vlim = abs(float(vlim_abs_env))
        else:
            vlim = DEFAULT_ESP_VLIM_ABS
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to read ESP cube for scaling (%s): %s", esp_cube_path, exc)
        vlim = DEFAULT_ESP_VLIM_ABS

    output_png.parent.mkdir(parents=True, exist_ok=True)
    code = r"""
import sys
xyz, esp_cube, out_png, transp_str, vlim_str, v1, v2, v3, v4, v5, v6 = sys.argv[1:12]
transp = float(transp_str)
vlim = float(vlim_str)
try:
    import pymol
    from pymol import cmd
except Exception as exc:
    print(f"PyMOL import failed: {exc}", file=sys.stderr)
    sys.exit(2)

def _style_scene():
    cmd.hide("everything", "mol")
    # Ball-and-stick like the dipole plot
    cmd.show("sticks", "mol")
    cmd.show("spheres", "mol")
    cmd.set("sphere_scale", 0.25)
    cmd.set("stick_radius", 0.15)
    cmd.set("stick_quality", 15)
    cmd.set("sphere_quality", 3)
    cmd.set("valence", 1)

    cmd.show("surface", "mol")
    cmd.set("transparency", transp, "mol")
    cmd.set("surface_quality", 1)

    cmd.color("gray50", "elem C")
    cmd.color("white", "elem H")
    cmd.color("slate", "elem N")
    cmd.color("red", "elem O")
    cmd.color("yellow", "elem S")

    cmd.bg_color("white")
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("ambient", 0.35)
    cmd.set("direct", 0.65)
    cmd.set("specular", 0.2)
    cmd.set("shininess", 10)
    cmd.set("depth_cue", 0)

def _apply_map():
    cmd.ramp_new("esp_ramp", "esp", [-vlim, 0.0, vlim], ["red", "white", "blue"])
    cmd.set("surface_color", "esp_ramp", "mol")
    # Hide the ramp object itself; otherwise PyMOL draws a horizontal legend bar in the render.
    try:
        cmd.disable("esp_ramp")
    except Exception:
        pass
    try:
        cmd.hide("everything", "esp_ramp")
    except Exception:
        pass

try:
    pymol.finish_launching(["pymol", "-cq"])
    cmd.reinitialize()

    cmd.load(xyz, "mol")
    _style_scene()

    cmd.load(esp_cube, "esp")
    _apply_map()

    cmd.orient("mol")
    cmd.zoom("all", buffer=2.5)

    # Render 6 views (3x2 montage assembled outside PyMOL)
    y_angles = [-60, 0, 60]
    x_angles = [15, -30]
    outs = [v1, v2, v3, v4, v5, v6]
    views = []
    k = 0
    for xa in x_angles:
        for ya in y_angles:
            views.append((outs[k], [("y", ya), ("x", xa)]))
            k += 1

    for out, turns in views:
        cmd.orient("mol")
        for axis, ang in turns:
            cmd.turn(axis, ang)
        cmd.zoom("all", buffer=2.5)
        cmd.ray(1200, 1000)
        cmd.png(out, dpi=300)

    # Leave a copy for compatibility; will be overwritten by montage afterwards.
    cmd.png(out_png, dpi=300)
except Exception as exc:
    print(f"PyMOL render failed: {exc}", file=sys.stderr)
    sys.exit(3)
finally:
    try:
        cmd.delete("all")
        cmd.reinitialize()
    except Exception:
        pass
"""
    proc = subprocess.run(
        [
            sys.executable,
            "-c",
            code,
            str(xyz_path),
            str(esp_cube_path),
            str(output_png),
            str(surface_transparency),
            str(vlim),
            str(output_png.with_name(output_png.stem + "_v1.png")),
            str(output_png.with_name(output_png.stem + "_v2.png")),
            str(output_png.with_name(output_png.stem + "_v3.png")),
            str(output_png.with_name(output_png.stem + "_v4.png")),
            str(output_png.with_name(output_png.stem + "_v5.png")),
            str(output_png.with_name(output_png.stem + "_v6.png")),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    v_paths = [
        output_png.with_name(output_png.stem + "_v1.png"),
        output_png.with_name(output_png.stem + "_v2.png"),
        output_png.with_name(output_png.stem + "_v3.png"),
        output_png.with_name(output_png.stem + "_v4.png"),
        output_png.with_name(output_png.stem + "_v5.png"),
        output_png.with_name(output_png.stem + "_v6.png"),
    ]
    if proc.returncode != 0 or not all(p.exists() for p in v_paths):
        logger.warning(
            "PyMOL surface render failed (rc=%s).\nstdout:\n%s\nstderr:\n%s",
            proc.returncode,
            (proc.stdout or "").strip(),
            (proc.stderr or "").strip(),
        )
        return None

    montage_png = output_png.with_name(output_png.stem + "_montage.png")
    montage = _make_grid_png(v_paths, montage_png, cols=3, pad_px=10)
    for p in v_paths:
        try:
            p.unlink()
        except Exception:
            pass
    if montage and montage.exists():
        try:
            montage.replace(output_png)
        except Exception:
            # fallback: overwrite by save
            try:
                from PIL import Image
                Image.open(montage).save(output_png)
            except Exception:
                pass
        try:
            montage_png.unlink()
        except Exception:
            pass

    import os

    if str(os.getenv("DELFIN_ESP_ADD_COLORBAR", "1")).lower() in {"0", "false", "no"}:
        return output_png

    colorbar_png = output_png.with_name(output_png.stem + "_colorbar.png")
    cb = _create_colorbar_png(colorbar_png, vlim=vlim)
    if cb and cb.exists():
        merged = _merge_with_colorbar(output_png, cb, output_png)
        try:
            colorbar_png.unlink()
        except Exception:
            pass
        return merged

    return output_png


def create_esp_slice_png(
    cube_path: Path,
    output_png: Path,
    slice_axis: str = "z",
) -> Optional[Path]:
    """Create a 2D mid-plane slice of an ESP cube with a colorbar."""
    try:
        import numpy as np
        import matplotlib.pyplot as plt
    except Exception as exc:  # noqa: BLE001
        logger.warning("matplotlib/numpy not available; cannot plot ESP: %s", exc)
        return None

    try:
        header, data = _read_cube(cube_path)
    except Exception as exc:  # noqa: BLE001
        logger.error("Failed to read cube %s: %s", cube_path, exc)
        return None

    axis = slice_axis.lower().strip()
    if axis not in {"x", "y", "z"}:
        axis = "z"

    if axis == "z":
        idx = header.nz // 2
        slice2d = data[idx, :, :]
        xlabel, ylabel = "X (Å)", "Y (Å)"
        x_step = (header.vx_bohr[0] ** 2 + header.vx_bohr[1] ** 2 + header.vx_bohr[2] ** 2) ** 0.5
        y_step = (header.vy_bohr[0] ** 2 + header.vy_bohr[1] ** 2 + header.vy_bohr[2] ** 2) ** 0.5
        extent = [0, header.nx * x_step * BOHR_TO_ANGSTROM, 0, header.ny * y_step * BOHR_TO_ANGSTROM]
        title = f"ESP slice ({cube_path.name}), z index {idx}"
    elif axis == "y":
        idx = header.ny // 2
        slice2d = data[:, idx, :]
        xlabel, ylabel = "X (Å)", "Z (Å)"
        x_step = (header.vx_bohr[0] ** 2 + header.vx_bohr[1] ** 2 + header.vx_bohr[2] ** 2) ** 0.5
        z_step = (header.vz_bohr[0] ** 2 + header.vz_bohr[1] ** 2 + header.vz_bohr[2] ** 2) ** 0.5
        extent = [0, header.nx * x_step * BOHR_TO_ANGSTROM, 0, header.nz * z_step * BOHR_TO_ANGSTROM]
        title = f"ESP slice ({cube_path.name}), y index {idx}"
    else:  # axis == "x"
        idx = header.nx // 2
        slice2d = data[:, :, idx]
        xlabel, ylabel = "Y (Å)", "Z (Å)"
        y_step = (header.vy_bohr[0] ** 2 + header.vy_bohr[1] ** 2 + header.vy_bohr[2] ** 2) ** 0.5
        z_step = (header.vz_bohr[0] ** 2 + header.vz_bohr[1] ** 2 + header.vz_bohr[2] ** 2) ** 0.5
        extent = [0, header.ny * y_step * BOHR_TO_ANGSTROM, 0, header.nz * z_step * BOHR_TO_ANGSTROM]
        title = f"ESP slice ({cube_path.name}), x index {idx}"

    finite = np.isfinite(slice2d)
    if not finite.any():
        logger.error("ESP slice contains no finite values: %s", cube_path)
        return None

    abs_p = np.abs(slice2d[finite])
    # Robust scaling; keep symmetric range around 0 for ESP.
    vlim = float(np.nanpercentile(abs_p, 99.0)) if abs_p.size else float(np.nanmax(abs_p))
    if not vlim or vlim <= 0:
        vlim = float(np.nanmax(abs_p)) if abs_p.size else 1.0

    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    im = ax.imshow(
        slice2d,
        origin="lower",
        cmap="coolwarm",
        vmin=-vlim,
        vmax=vlim,
        extent=extent,
        interpolation="nearest",
        aspect="equal",
    )
    # Remove axes/labels to avoid matplotlib "bars" under the plot in the final report image.
    ax.set_axis_off()

    import os
    if str(os.getenv("DELFIN_ESP_ADD_COLORBAR", "1")).lower() not in {"0", "false", "no"}:
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        try:
            cbar.outline.set_visible(False)
        except Exception:
            pass
        for spine in cbar.ax.spines.values():
            spine.set_visible(False)
        cbar.ax.tick_params(labelsize=8, length=0)
        cbar.set_label("a.u.", fontsize=8, labelpad=2)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    fig.savefig(str(output_png), dpi=300, bbox_inches="tight", pad_inches=0.0)
    plt.close(fig)
    logger.info("ESP slice PNG saved to %s", output_png)
    return output_png if output_png.exists() else None


def generate_esp_png(workspace_root: Path) -> Optional[Path]:
    """Generate an ESP visualization PNG for S0 (3D isosurface if possible, else 2D slice)."""
    return generate_esp_png_for_state(workspace_root, "S0")


def _state_file_candidates(state_name: str, ext: str) -> list[str]:
    state = state_name.upper()
    candidates = [
        f"{state}.{ext}",
        f"{state}_second_deltaSCF.{ext}",
        f"{state}_TDDFT.{ext}",
        f"{state}_opt.{ext}",
    ]
    if state == "S0":
        candidates.append(f"initial.{ext}")
    return candidates


def generate_esp_png_for_state(workspace_root: Path, state_name: str) -> Optional[Path]:
    """Generate an ESP visualization PNG for a specific state (3D isosurface if possible, else 2D slice)."""
    gbw = _find_first_existing(workspace_root, _state_file_candidates(state_name, "gbw"))
    if not gbw:
        logger.warning("No GBW file found for %s; skipping ESP plot", state_name)
        return None

    # Locate / generate ESP cube (in GBW directory)
    start_esp = time.time()
    esp_cube = _pick_esp_cube(gbw.parent)
    if not esp_cube:
        ok = _orca_plot_generate_esp_cube(gbw)
        if ok:
            esp_cube = _pick_esp_cube(gbw.parent, modified_after=start_esp - 1.0) or _pick_esp_cube(gbw.parent)
    if not esp_cube:
        logger.warning("No ESP cube found/generated next to %s", gbw)
        return None

    # Prefer 3D render if we have an XYZ; color the molecular surface by the ESP map.
    xyz = _find_first_existing(workspace_root, _state_file_candidates(state_name, "xyz"))
    if xyz:
        out_3d = workspace_root / f"Electrostatic_Potential_{state_name.upper()}.png"
        rendered = create_esp_molecular_surface_png(xyz, esp_cube, out_3d)
        if rendered:
            return rendered

    # Fallback: 2D slice with colorbar
    output_png = workspace_root / f"Electrostatic_Potential_{state_name.upper()}.png"
    return create_esp_slice_png(esp_cube, output_png, slice_axis="z")
