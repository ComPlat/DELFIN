from __future__ import annotations

from html import escape


def controls_html(bio: dict) -> str:
    options = '<option value="__all__">All chains</option>' + ''.join(
        f'<option value="{escape(chain, quote=True)}">{escape(chain or "(blank)")}</option>'
        for chain in bio["chains"]
    )
    return (
        '<span class="delfin-bio-controls" style="display:flex;gap:8px;align-items:center;flex-wrap:wrap;">'
        f'<label>Chain <select class="delfin-bio-chain">{options}</select></label>'
        '<label><input type="checkbox" class="delfin-bio-waters"> Waters</label></span>'
    )
