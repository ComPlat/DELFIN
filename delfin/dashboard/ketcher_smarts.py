"""What was drawn, as a SMARTS that RDKit means the same thing by.

Ketcher can write a SMARTS itself -- ``ketcher.getSmarts()``, which is Indigo
writing it -- and for the left of the arrow that is enough: Indigo's
``[#6;a;H1:1]`` matches exactly what a hand-written ``[cH:1]`` matches.

The right of the arrow is where it stops being enough.  A product template is
not matched, it is *built*, and RDKit builds it from the atom it is given.  A
query atom carries no element aromaticity and no hydrogen count, so Indigo's

    [#6;a;H1:1]:...:[#6;a:6] >> [#8;a:2]1[#6;a:3][#6;a:4][#6;a:5][#6;a:6]1

turns benzene into *tetrahydrofuran*, where the same rule written by hand,

    [c&H1:1]:...:[c:6] >> [o:2]1[c:3][c:4][c:5][c:6]1

turns it into furan.  Nothing errors; the wrong molecule is simply what comes
out.  So the product side is rebuilt here out of real atoms -- element,
aromaticity, charge and map -- before the reaction is assembled.

The one thing deliberately *not* carried over is the hydrogen count.  Every
mapped atom is written in brackets, and a bracket atom in SMILES states its
hydrogens; let RDKit work them out and ``[c:3]`` becomes ``[cH:3]``, which
pins the product at exactly one hydrogen and quietly stops the rule from
firing on a fused position.  ``SetNoImplicit`` before the write is what keeps
it unspecified.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple

from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions, rdFMCS

__all__ = [
    'split_reaction', 'concretize_product', 'complete_atom_maps',
    'prettify', 'normalize_reaction_smarts', 'normalize_query_smarts',
    'describe', 'inspect', 'reaction_smarts_from_rxn_block',
    'query_smarts_from_molblock', 'rxn_block_from_smarts', 'widen_kekule',
    'merge_maps',
    'absorb_hydrogens',
    'trial_on_seed',
    'survives_the_editor',
]

RDLogger.DisableLog('rdApp.*')

#: Aromatic and aliphatic shorthands, so what comes back reads like what a
#: chemist would have typed.  Applied only when the round trip proves the
#: shorter form means the same thing, so a wrong guess here costs nothing.
_SHORTHAND = [
    ('#6&a', 'c'), ('#7&a', 'n'), ('#8&a', 'o'),
    ('#16&a', 's'), ('#15&a', 'p'), ('#5&a', 'b'),
    ('#6&A', 'C'), ('#7&A', 'N'), ('#8&A', 'O'),
    ('#16&A', 'S'), ('#15&A', 'P'),
]

_BARE = re.compile(r'\[([cnops])\]')


def _depth_split(text: str, sep: str = '>') -> List[str]:
    """Split on *sep*, but not inside brackets.

    ``>`` is a field separator in a reaction SMARTS and also a character that
    turns up inside a recursive query, so counting depth is the difference
    between reading ``A>>B`` and cutting ``[$(C>...)]`` in half.
    """
    out: List[str] = []
    depth = 0
    start = 0
    for i, ch in enumerate(text):
        if ch in '[(':
            depth += 1
        elif ch in '])':
            depth -= 1
        elif ch == sep and depth <= 0:
            out.append(text[start:i])
            start = i + 1
    out.append(text[start:])
    return out


def _atom_blocks(text: str):
    """Walk the bracket atoms, yielding ``(start, end)`` over the inside.

    Depth-aware rather than a regular expression: a recursive query carries
    brackets of its own, and ``[$([#6])]`` cut at the first ``]`` is two
    halves of nothing.
    """
    depth = 0
    start = -1
    for i, ch in enumerate(text):
        if ch == '[':
            if depth == 0:
                start = i + 1
            depth += 1
        elif ch == ']':
            depth -= 1
            if depth == 0 and start >= 0:
                yield start, i
                start = -1


_MAP_INSIDE = re.compile(r':(\d+)')


def _maps_in(inside: str):
    """The map numbers of one atom, ignoring anything inside ``$(...)``.

    A recursive query holds bonds, and a bond is written ``:`` too.  Reading
    ``[#6:2;$([*,#1]=,#,:[*,#1])]`` without counting depth finds the aromatic
    bond in the recursion and gives up on the atom -- which is how Ketcher's
    Unsaturated came through as something RDKit could not read at all.
    """
    found, depth = [], 0
    for hit in re.finditer(r'[()]|:(\d+)', inside):
        token = hit.group(0)
        if token == '(':
            depth += 1
        elif token == ')':
            depth -= 1
        elif depth == 0:
            found.append(hit)
    return found


def repair_atom_maps(text: str) -> str:
    """Move an atom map back to the end of its atom, where SMARTS wants it.

    Ketcher writes ``[#6:3;v2]``.  The map has to come last -- ``[#6;v2:3]``
    -- and RDKit is right to refuse the other order, so a drawing that picked
    up a valence query on the way through the editor arrives unreadable.  It
    is the same atom either way, so it is put right rather than reported.

    It also repeats the map on every term of a list rather than writing it
    once.  A halogen drawn as Ketcher's X comes back as
    ``[F:2,Cl:2,Br:2,I:2,At:2]`` and a metal drawn as M as a seventeen-term
    negation carrying ``:2`` seventeen times -- neither parses, so X and M
    reached nothing at all before this.  Repeats of one number collapse to a
    single map at the end; two *different* numbers in one atom are a genuine
    ambiguity and are left alone to be reported.

    A ``:`` inside a recursive query is a bond, not a map, so the search for
    one counts bracket depth rather than skipping such atoms -- Ketcher writes
    its Unsaturated as a recursive query with the map in the wrong place, and
    skipping it left the whole rule unreadable.
    """
    original = str(text or '')
    out = []
    last = 0
    for start, end in _atom_blocks(original):
        inside = original[start:end]
        found = _maps_in(inside)
        if not found:
            continue
        numbers = {hit.group(1) for hit in found}
        if len(numbers) != 1:
            continue
        if len(found) == 1 and found[0].end() == len(inside):
            continue
        moved, cut = [], 0
        for hit in found:
            moved.append(inside[cut:hit.start()])
            cut = hit.end()
        moved.append(inside[cut:])
        rebuilt = ''.join(moved).rstrip(';&,')
        out.append(original[last:start])
        out.append(f'{rebuilt}:{numbers.pop()}')
        last = end
    if not out:
        return original
    out.append(original[last:])
    return ''.join(out)


def split_reaction(raw: str) -> Optional[Tuple[str, str, str]]:
    """``(reactants, agents, products)``, whichever way the arrow was written.

    Ketcher writes a plain two-field reaction; a three-field one arrives when
    something was drawn over the arrow.  Both are accepted, and the agents are
    kept so they can be reported rather than silently dropped.
    """
    text = ' '.join(str(raw or '').split())
    if not text:
        return None
    fields = _depth_split(text, '>')
    if len(fields) == 3:
        return fields[0].strip(), fields[1].strip(), fields[2].strip()
    if len(fields) == 2:
        return fields[0].strip(), '', fields[1].strip()
    return None


def _is_aromatic(atom) -> bool:
    """Whether this atom is aromatic, however it came to be written.

    A query atom's own flag is never set -- ``[#8;a]`` parses to an AtomAnd
    over an atomic number and an aromaticity test, and the atom itself stays
    unaromatic -- so the query has to be read rather than the atom.
    """
    if atom.GetIsAromatic():
        return True
    if not atom.HasQuery():
        return False
    try:
        return 'AtomIsAromatic 1 = val' in atom.DescribeQuery()
    except Exception:                                       # noqa: BLE001
        return False


def concretize_product(template: Any) -> Optional[str]:
    """Rewrite a product template out of real atoms.

    Takes a SMARTS string or a mol; gives back the template as text RDKit will
    build the intended molecule from.  See the module docstring for why this
    exists at all.
    """
    query = (Chem.MolFromSmarts(template)
             if isinstance(template, str) else template)
    if query is None:
        return None

    built = Chem.RWMol()
    seats: Dict[int, int] = {}
    for atom in query.GetAtoms():
        fresh = Chem.Atom(atom.GetAtomicNum())
        fresh.SetFormalCharge(atom.GetFormalCharge())
        fresh.SetAtomMapNum(atom.GetAtomMapNum())
        fresh.SetIsAromatic(_is_aromatic(atom))
        # Left unspecified on purpose -- the module docstring says why.
        fresh.SetNoImplicit(True)
        fresh.SetNumExplicitHs(0)
        seats[atom.GetIdx()] = built.AddAtom(fresh)

    for bond in query.GetBonds():
        here, there = seats[bond.GetBeginAtomIdx()], seats[bond.GetEndAtomIdx()]
        both_aromatic = (built.GetAtomWithIdx(here).GetIsAromatic()
                         and built.GetAtomWithIdx(there).GetIsAromatic())
        kind = bond.GetBondType()
        try:
            described = bond.DescribeQuery()
        except Exception:                                   # noqa: BLE001
            described = ''
        # A ring closure Indigo wrote as single-or-aromatic, between two atoms
        # it also called aromatic, is an aromatic bond.  Left as SINGLE it is
        # a saturated ring with aromatic atoms in it, which sanitizes into
        # something nobody drew.
        if both_aromatic and ('SingleOrAromatic' in described
                              or kind in (Chem.BondType.SINGLE,
                                          Chem.BondType.UNSPECIFIED,
                                          Chem.BondType.AROMATIC)):
            kind = Chem.BondType.AROMATIC
        built.AddBond(here, there, kind)
        if kind == Chem.BondType.AROMATIC:
            built.GetBondBetweenAtoms(here, there).SetIsAromatic(True)

    # An atom in an aromatic bond is an aromatic atom, whether or not whatever
    # drew it said so.  A structure that came the RXN way rather than the
    # SMARTS way carries its aromaticity only on the bonds -- the atoms arrive
    # as bare ``[#6]`` -- and without this the flag never reaches the product.
    for atom in built.GetAtoms():
        if atom.GetIsAromatic():
            continue
        if any(b.GetBondType() == Chem.BondType.AROMATIC
               for b in atom.GetBonds()):
            atom.SetIsAromatic(True)

    out = built.GetMol()
    try:
        out.UpdatePropertyCache(strict=False)
    except Exception:                                       # noqa: BLE001
        pass
    try:
        return Chem.MolToSmiles(out, canonical=False)
    except Exception:                                       # noqa: BLE001
        return None


def _plain(query) -> Optional[Any]:
    """A real molecule with the query's atom order, for the MCS to chew on.

    ``FindMCS`` wants molecules, not patterns, and the indices have to line up
    with the pattern's so a match can be carried back as a map number.
    """
    try:
        built = Chem.RWMol()
        for atom in query.GetAtoms():
            fresh = Chem.Atom(atom.GetAtomicNum())
            fresh.SetIsAromatic(_is_aromatic(atom))
            fresh.SetNoImplicit(True)
            built.AddAtom(fresh)
        for bond in query.GetBonds():
            kind = bond.GetBondType()
            if kind == Chem.BondType.UNSPECIFIED:
                kind = Chem.BondType.SINGLE
            built.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), kind)
        out = built.GetMol()
        out.UpdatePropertyCache(strict=False)
        return out
    except Exception:                                       # noqa: BLE001
        return None


def complete_atom_maps(reactant, product) -> List[int]:
    """Fill in the maps nobody set, and leave the ones somebody did.

    Ketcher's mapping tool is fiddly and a forgotten map does not fail -- it
    changes the reaction into a different one that still runs.  So the common
    core of the two sides is found and the pairs in it that are unmapped on
    *both* sides get the next free numbers.  A map set by hand always wins:
    this only ever fills holes, so a deliberate deletion (mapped left, absent
    right) is never papered over.

    Returns the numbers it handed out, so the panel can say what it did.
    """
    given: List[int] = []
    used = {a.GetAtomMapNum() for a in reactant.GetAtoms() if a.GetAtomMapNum()}
    used |= {a.GetAtomMapNum() for a in product.GetAtoms() if a.GetAtomMapNum()}

    free_left = [a for a in reactant.GetAtoms() if not a.GetAtomMapNum()]
    free_right = [a for a in product.GetAtoms() if not a.GetAtomMapNum()]
    if not free_left or not free_right:
        return given

    here, there = _plain(reactant), _plain(product)
    if here is None or there is None:
        return given

    def pair(atoms_match):
        """One MCS pass; returns the index pairs it can line up."""
        try:
            found = rdFMCS.FindMCS(
                [here, there],
                atomCompare=atoms_match,
                # The point of the drawing is that the bonds change; comparing
                # them would refuse the very cores this is meant to line up.
                bondCompare=rdFMCS.BondCompare.CompareAny,
                ringMatchesRingOnly=False, completeRingsOnly=False,
                timeout=5,
            )
            if found.canceled or not found.smartsString:
                return []
            core = Chem.MolFromSmarts(found.smartsString)
            if core is None:
                return []
            left_hit = here.GetSubstructMatch(core)
            right_hit = there.GetSubstructMatch(core)
        except Exception:                                   # noqa: BLE001
            return []
        if not left_hit or not right_hit or len(left_hit) != len(right_hit):
            return []
        return list(zip(left_hit, right_hit))

    nxt = 1

    def hand_out(pairs):
        nonlocal nxt
        for l_idx, r_idx in pairs:
            l_atom = reactant.GetAtomWithIdx(l_idx)
            r_atom = product.GetAtomWithIdx(r_idx)
            if l_atom.GetAtomMapNum() or r_atom.GetAtomMapNum():
                continue
            while nxt in used:
                nxt += 1
            l_atom.SetAtomMapNum(nxt)
            r_atom.SetAtomMapNum(nxt)
            used.add(nxt)
            given.append(nxt)

    # Elements first, which is the safe pairing.
    hand_out(pair(rdFMCS.AtomCompare.CompareElements))
    # Then, for whatever is still unmapped on both sides, elements are allowed
    # to differ.  A morphing rule exists to change an element, so the atom the
    # map matters most for -- the carbon that becomes a nitrogen -- is exactly
    # the one an element-comparing MCS will not pair.  Unmapped it would be
    # deleted and its replacement built loose, which is a rule that matches and
    # makes nothing.
    if (any(not a.GetAtomMapNum() for a in reactant.GetAtoms())
            and any(not a.GetAtomMapNum() for a in product.GetAtoms())):
        hand_out(pair(rdFMCS.AtomCompare.CompareAny))
    return given


#: What a rewritten pattern is tried against to see whether it still means the
#: same thing.  Small, aromatic and aliphatic, with the heteroatoms the rules
#: in this tab actually swap in and out.
_PROBES = (
    'c1ccccc1', 'Cc1ccccc1', 'c1ccc2ccccc2c1', 'c1ccncc1', 'c1ccoc1',
    'c1ccsc1', 'Oc1ccccc1', 'C1CCCCC1', 'CCCCCC', 'CC(=O)O', 'C=Cc1ccccc1',
)
_PROBE_MOLS: List[Any] = []


#: A bond that is either the Kekule form or the aromatic one.  Taken off a
#: parsed template because a bond query cannot be built from Python directly --
#: ``QueryBond`` exposes no ``GetQuery`` -- and copied on with ``ReplaceBond``.
_SINGLE_OR_AROMATIC = Chem.MolFromSmarts('[#6]-,:[#6]').GetBondWithIdx(0)
_DOUBLE_OR_AROMATIC = Chem.MolFromSmarts('[#6]=,:[#6]').GetBondWithIdx(0)


#: One aromatic atom of each element a drawn fragment is likely to sit on, to
#: ask a query atom whether it would accept an aromatic partner at all.
_AROMATIC_SOURCES = ('c1ccccc1', 'c1ccncc1', 'c1ccoc1', 'c1ccsc1', 'c1cc[se]c1')
_AROMATIC_ATOMS: List[Any] = []


def _aromatic_probes() -> List[Any]:
    if not _AROMATIC_ATOMS:
        for smiles in _AROMATIC_SOURCES:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            for atom in mol.GetAtoms():
                if atom.GetIsAromatic():
                    _AROMATIC_ATOMS.append(atom)
    return _AROMATIC_ATOMS


def _may_be_aromatic(atom) -> bool:
    """Whether this query atom would accept an aromatic atom at all.

    ``[#6]`` would, ``[C]`` and ``[#6;A]`` would not -- the difference between
    a carbon drawn without saying which kind and one drawn as explicitly
    aliphatic, which is exactly the line a widening must not cross.
    """
    if not atom.HasQuery():
        return bool(atom.GetIsAromatic())
    for probe in _aromatic_probes():
        try:
            if atom.Match(probe):
                return True
        except Exception:                                   # noqa: BLE001
            return False
    return False


def absorb_hydrogens(smarts: str) -> Tuple[str, int]:
    """Turn a hydrogen drawn as an atom into the hydrogen count it was meant as.

    Putting an H on a free position is the natural thing to draw, and in a
    query it is the one thing that cannot work: an explicit ``[#1]`` matches an
    explicit hydrogen, and a molecule RDKit has read carries its hydrogens
    implicitly.  Anthracene has none at all to match, so a rule drawn that way
    finds nothing -- measured: 0 matches with the H drawn, 8 with it stated as
    ``H1`` on the atom it hangs off.

    Each drawn hydrogen is removed and counted onto its neighbour instead.  The
    neighbour is rebuilt from its own SMARTS text with ``;H<n>`` spliced in
    ahead of the map, so whatever else was set on it survives.
    """
    text = str(smarts or '')
    query = Chem.MolFromSmarts(text)
    if query is None:
        return text, 0
    loose: Dict[int, int] = {}
    for atom in query.GetAtoms():
        if atom.GetAtomicNum() != 1 or atom.GetDegree() != 1:
            continue
        host = atom.GetNeighbors()[0]
        if host.GetAtomicNum() == 1:
            continue
        loose[host.GetIdx()] = loose.get(host.GetIdx(), 0) + 1
    if not loose:
        return text, 0

    out = Chem.RWMol(query)
    for index, count in loose.items():
        written = query.GetAtomWithIdx(index).GetSmarts()
        if not (written.startswith('[') and written.endswith(']')):
            written = f'[{written}]'
        inside = written[1:-1]
        head, mark, tail = inside.rpartition(':')
        if mark and tail.isdigit():
            fresh = f'[{head};H{count}:{tail}]'
        else:
            fresh = f'[{inside};H{count}]'
        made = Chem.MolFromSmarts(fresh)
        if made is None or made.GetNumAtoms() != 1:
            return text, 0
        try:
            out.ReplaceAtom(index, made.GetAtomWithIdx(0), False, True)
        except Exception:                                   # noqa: BLE001
            return text, 0
    for index in sorted((a.GetIdx() for a in query.GetAtoms()
                         if a.GetAtomicNum() == 1 and a.GetDegree() == 1),
                        reverse=True):
        try:
            out.RemoveAtom(index)
        except Exception:                                   # noqa: BLE001
            return text, 0
    try:
        return Chem.MolToSmarts(out.GetMol()), sum(loose.values())
    except Exception:                                       # noqa: BLE001
        return text, 0


def widen_kekule(smarts: str, *, everywhere: bool = True) -> Tuple[str, int]:
    """Let a ring drawn in Kekule form match the aromatic ring it stands for.

    Ketcher draws benzene the way a chemist does, with three alternating
    double bonds, and Indigo writes that out as it was drawn:
    ``[#6:1]1-[#6:2]=[#6:3]-...``.  The seed on the other side of the match is
    parsed by RDKit and comes out *aromatic* -- every bond ``AROMATIC``, no
    single or double anywhere -- and a ``-`` query does not match an aromatic
    bond, nor a ``=``.  So the rule matched nothing at all and the tab said
    "Reaction SMARTS produced no products", with a drawing that was right.

    Each ring bond that aromaticity perception calls aromatic is widened to
    ``-,:`` or ``=,:``, which matches the ring in either form.  Only those:
    cyclohexa-1,4-diene is drawn with the same two bond types and perceives as
    nothing, so it stays exactly as narrow as it was drawn.

    A ring is the case this can decide on its own.  A fragment *cut out* of a
    ring cannot be: ``[#6:1](~[#6:2])=[#6]~[#6]`` is four atoms of a benzene
    with one bond drawn double, and nothing in it says the double bond stands
    for an aromatic one.  So with *everywhere* -- which the panel offers as a
    checkbox, on by default because ChemDarwin's seeds are aromatic scaffolds
    -- every drawn single or double bond is widened as well, but only between
    two atoms that would accept an aromatic partner.  A bond drawn between
    explicitly aliphatic atoms is left alone either way.

    Returns the SMARTS and how many bonds were widened, so the panel can say
    that it happened rather than quietly changing what somebody drew.
    """
    text = str(smarts or '')
    query = Chem.MolFromSmarts(text)
    if query is None:
        return text, 0
    plain = _plain(query)
    if plain is None:
        return text, 0
    try:
        Chem.SanitizeMol(plain)
    except Exception:                                       # noqa: BLE001
        # No perception, so no ring to decide on -- but the *everywhere* pass
        # asks the atoms, not the ring, and still has something to say.
        if not everywhere:
            return text, 0

    out = Chem.RWMol(query)
    widened = 0
    for index, reference in enumerate(plain.GetBonds()):
        drawn = query.GetBondWithIdx(index)
        kind = drawn.GetBondType()
        if kind not in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE):
            continue
        if not reference.GetIsAromatic():
            if not everywhere:
                continue
            if not (_may_be_aromatic(drawn.GetBeginAtom())
                    and _may_be_aromatic(drawn.GetEndAtom())):
                continue
        template = (_DOUBLE_OR_AROMATIC if kind == Chem.BondType.DOUBLE
                    else _SINGLE_OR_AROMATIC)
        try:
            out.ReplaceBond(index, template, True)
        except Exception:                                   # noqa: BLE001
            return text, 0
        widened += 1
    if not widened:
        return text, 0
    try:
        return Chem.MolToSmarts(out.GetMol()), widened
    except Exception:                                       # noqa: BLE001
        return text, 0


def _probes() -> List[Any]:
    if not _PROBE_MOLS:
        for smiles in _PROBES:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                _PROBE_MOLS.append(mol)
    return _PROBE_MOLS


def _behaviour(text: str, reaction: bool) -> Optional[frozenset]:
    """What *text* does to the probes -- None if it cannot be read at all.

    Comparing what RDKit *writes* would be no use here: ``[#6&a]`` and ``[c]``
    are the same atom and canonicalise to different strings, which is the
    whole reason the short spelling is worth having.  So the two are compared
    by what they match, or by what they make.
    """
    try:
        if reaction:
            built = rdChemReactions.ReactionFromSmarts(text)
            if built is None:
                return None
            made = set()
            for probe in _probes():
                for group in built.RunReactants((probe,)):
                    for mol in group:
                        try:
                            copy = Chem.Mol(mol)
                            Chem.SanitizeMol(copy)
                            made.add(Chem.MolToSmiles(copy))
                        except Exception:               # noqa: BLE001
                            made.add('?')
            return frozenset(made)
        query = Chem.MolFromSmarts(text)
        if query is None:
            return None
        hits = set()
        for i, probe in enumerate(_probes()):
            for hit in probe.GetSubstructMatches(query, uniquify=True):
                hits.add((i, hit))
        return frozenset(hits)
    except Exception:                                       # noqa: BLE001
        return None


def prettify(text: str, *, reaction: bool = True) -> str:
    """``[#6&a&H1:1]`` back to ``[c&H1:1]``, when it is safe to do so.

    Purely how it is spelled.  RDKit writes every element as an atomic number,
    and the rules already sitting in the box beside this one were written by
    hand; a box that mixes the two dialects is a box nobody can read.

    The shortened form is kept only if it still matches and makes exactly what
    the long one did, measured against `_PROBES` -- so a substitution that
    landed somewhere unintended is dropped rather than shipped.
    """
    original = str(text or '')
    if not original:
        return original
    shorter = original
    for verbose, brief in _SHORTHAND:
        shorter = shorter.replace(verbose, brief)
    shorter = _BARE.sub(r'\1', shorter)
    if shorter == original:
        return original
    was = _behaviour(original, reaction)
    now = _behaviour(shorter, reaction)
    if was is not None and was == now:
        return shorter
    return original


def _maps(mol) -> List[int]:
    return sorted({a.GetAtomMapNum() for a in mol.GetAtoms()
                   if a.GetAtomMapNum()})


def _listed(numbers: List[int]) -> str:
    return ', '.join(str(n) for n in numbers)


def describe(outcome: Dict[str, Any]) -> str:
    """The one line under the preview: what this rule will actually do."""
    if not outcome.get('ok'):
        return str(outcome.get('status') or 'Cannot be read.')
    bits = [f"{outcome['reactants']} reactant"
            f"{'s' if outcome['reactants'] != 1 else ''}",
            f"{outcome['products']} product"
            f"{'s' if outcome['products'] != 1 else ''}"]
    dropped = outcome.get('deleted_maps') or []
    if dropped:
        bits.append(f"map {_listed(dropped)} is deleted")
    if outcome.get('deleted_unmapped'):
        bits.append(f"{outcome['deleted_unmapped']} unmapped atom"
                    f"{'s' if outcome['deleted_unmapped'] != 1 else ''} on the "
                    f"left, deleted")
    if outcome.get('new_atoms'):
        bits.append(f"{outcome['new_atoms']} new atom"
                    f"{'s' if outcome['new_atoms'] != 1 else ''}")
    if outcome.get('auto_maps'):
        if outcome.get('drew_maps') is False:
            bits.append('NOTHING was mapped in the drawing — the mapping below '
                        'is a guess; join the atoms that stay with the '
                        'Reaction Mapping Tool')
        else:
            bits.append(f"map {_listed(outcome['auto_maps'])} "
                        f"filled in automatically")
    if outcome.get('kekule'):
        bits.append(f"{outcome['kekule']} drawn bond"
                    f"{'s' if outcome['kekule'] != 1 else ''} also match "
                    f"aromatic")
    if outcome.get('drawn_h'):
        bits.append(f"{outcome['drawn_h']} drawn hydrogen"
                    f"{'s' if outcome['drawn_h'] != 1 else ''} read as an "
                    f"H count")
    if outcome.get('agents'):
        bits.append('what was drawn over the arrow is ignored')
    return ' · '.join(bits)


def normalize_reaction_smarts(raw: str, *,
                              aromatic: bool = True) -> Dict[str, Any]:
    """A drawn reaction, as a line the Reaction SMARTS box can hold.

    ``level`` is what the panel colours by, and it is not the same question as
    ``ok``.  ``rxn.Validate`` calls "mapped in the reactants, not in the
    products" an *error*, and that is exactly what the furan rule does on
    purpose -- atom 1 leaves.  A rule RDKit builds and runs is therefore never
    refused here; it is reported, in the words of what it will do.
    """
    empty: Dict[str, Any] = {'ok': False, 'smarts': '', 'level': 'bad',
                             'auto_maps': [], 'deleted_maps': []}
    parts = split_reaction(repair_atom_maps(raw))
    if parts is None:
        return dict(empty, status='No reaction arrow in the drawing -- put '
                                  'one there with the Reaction Arrow Tool.')
    left, agents, right = parts
    if not left or not right:
        return dict(empty, status='One side of the arrow is empty.')

    # A hydrogen drawn as an atom matches an explicit hydrogen, and the
    # molecules this rule will meet carry theirs implicitly.  Both sides, so
    # the two templates keep the same atoms.
    left, drawn_h = absorb_hydrogens(left)
    right, right_h = absorb_hydrogens(right)
    drawn_h += right_h
    # Ketcher draws benzene in Kekule form and Indigo writes it out that way,
    # while every molecule this rule will meet has been through RDKit's
    # aromaticity perception.  Without this the two never touch.
    left, kekule = widen_kekule(left, everywhere=aromatic)
    reactant = Chem.MolFromSmarts(left)
    if reactant is None:
        return dict(empty, status=f'The reactant side cannot be read: {left}')
    product = Chem.MolFromSmarts(right)
    if product is None:
        return dict(empty, status=f'The product side cannot be read: {right}')

    # Whether the *drawing* said anything about what stays, before the MCS is
    # let near it.  A rule mapped entirely by guesswork is a different claim
    # from one where the holes were filled, and it has to read differently.
    drew_maps = bool(_maps(reactant) or _maps(product))
    auto = complete_atom_maps(reactant, product)

    try:
        left_text = Chem.MolToSmarts(reactant)
    except Exception as exc:                                # noqa: BLE001
        return dict(empty, status=f'The reactant side could not be written: {exc}')
    right_text = concretize_product(product)
    if right_text is None:
        return dict(empty, status='The product side could not be written.')

    smarts = f'{left_text}>>{right_text}'
    built = None
    try:
        built = rdChemReactions.ReactionFromSmarts(smarts)
    except Exception:                                       # noqa: BLE001
        built = None
    if built is None:
        return dict(empty, status=f'RDKit cannot build a reaction from this: {smarts}')
    try:
        built.Initialize()
        built.Validate(silent=True)
    except Exception:                                       # noqa: BLE001
        pass

    left_maps, right_maps = _maps(reactant), _maps(product)
    outcome: Dict[str, Any] = {
        'ok': True,
        'smarts': prettify(smarts, reaction=True),
        'level': 'ok',
        'reactants': built.GetNumReactantTemplates(),
        'products': built.GetNumProductTemplates(),
        'auto_maps': auto,
        'deleted_maps': [n for n in left_maps if n not in right_maps],
        'deleted_unmapped': sum(1 for a in reactant.GetAtoms()
                                if not a.GetAtomMapNum()),
        'new_atoms': sum(1 for a in product.GetAtoms()
                         if not a.GetAtomMapNum()),
        'agents': agents,
        'kekule': kekule,
        'drawn_h': drawn_h,
        'drew_maps': drew_maps,
    }
    if not left_maps and not right_maps:
        outcome['level'] = 'note'
        outcome['status'] = ('Nothing is mapped, and there was nothing to '
                             'line up -- join the atoms that stay the same '
                             'with the Reaction Mapping Tool.')
        return outcome
    outcome['status'] = describe(outcome)
    if auto or outcome['agents'] or kekule or drawn_h:
        outcome['level'] = 'note'
    return outcome


def normalize_query_smarts(raw: str, *,
                           aromatic: bool = True) -> Dict[str, Any]:
    """A drawn fragment, as a line the Forbidden or Protected box can hold.

    No arrow: these are matched against a molecule, never applied to one.  Map
    numbers are dropped -- ``GetSubstructMatches`` ignores them anyway, and a
    stray ``:1`` in a filter line only ever reads as a half-finished reaction.
    """
    text = ' '.join(repair_atom_maps(raw).split())
    if not text:
        return {'ok': False, 'smarts': '', 'level': 'bad',
                'status': 'The canvas is empty.'}
    if len(_depth_split(text, '>')) > 1:
        return {'ok': False, 'smarts': '', 'level': 'bad',
                'status': ('A pattern belongs here, not a reaction -- take '
                           'the arrow out of the drawing.')}
    text, drawn_h = absorb_hydrogens(text)
    text, kekule = widen_kekule(text, everywhere=aromatic)
    query = Chem.MolFromSmarts(text)
    if query is None:
        return {'ok': False, 'smarts': '', 'level': 'bad',
                'status': f'RDKit cannot read this pattern: {text}'}
    for atom in query.GetAtoms():
        atom.SetAtomMapNum(0)
    try:
        written = Chem.MolToSmarts(query)
    except Exception as exc:                                # noqa: BLE001
        return {'ok': False, 'smarts': '', 'level': 'bad',
                'status': f'The pattern could not be written: {exc}'}
    count = query.GetNumAtoms()
    outcome = {'ok': True, 'smarts': prettify(written, reaction=False),
               'level': 'ok',
               'status': f"{count} atom{'s' if count != 1 else ''} · valid pattern"}
    if drawn_h:
        outcome['level'] = 'note'
        outcome['status'] += (f" · {drawn_h} drawn hydrogen"
                              f"{'s' if drawn_h != 1 else ''} read as an "
                              f"H count")
    if kekule:
        outcome['level'] = 'note'
        outcome['status'] += (f" · {kekule} drawn bond"
                              f"{'s' if kekule != 1 else ''} also match "
                              f"aromatic")
    # Two things drawn side by side are one pattern that wants both at once,
    # which is rarely what somebody drawing two of them meant.  The box does
    # have an "any of these", but it is the ``;`` between separate patterns,
    # not a second fragment inside one.
    try:
        pieces = len(Chem.GetMolFrags(query))
    except Exception:                                       # noqa: BLE001
        pieces = 1
    if pieces > 1:
        outcome['level'] = 'note'
        outcome['status'] = (
            f'{count} atoms in {pieces} separate fragments -- as one pattern '
            f'that requires all of them at once. For "any of these", draw them '
            f'one at a time; they are joined with ;.')
    return outcome


def inspect(text: str, *, reaction: bool = True) -> Dict[str, Any]:
    """Whether this line is usable, exactly as it stands.

    Nothing is rewritten.  The preview box beside the editor is editable, and
    a box that rewrites what is in it while somebody is typing cannot be typed
    in -- so reading the drawing normalises, and this only ever judges.
    """
    line = str(text or '').strip()
    if not line:
        return {'ok': False, 'level': 'bad', 'status': ''}
    if not reaction:
        query = Chem.MolFromSmarts(line)
        if query is None:
            return {'ok': False, 'level': 'bad',
                    'status': 'RDKit cannot read this pattern.'}
        count = query.GetNumAtoms()
        return {'ok': True, 'level': 'ok',
                'status': f"{count} atom{'s' if count != 1 else ''} · valid pattern"}

    if '>>' not in line:
        return {'ok': False, 'level': 'bad',
                'status': 'A rule needs the form A>>B.'}
    try:
        built = rdChemReactions.ReactionFromSmarts(line)
    except Exception:                                       # noqa: BLE001
        built = None
    if built is None:
        return {'ok': False, 'level': 'bad',
                'status': 'RDKit cannot build a reaction from this.'}
    try:
        built.Initialize()
        built.Validate(silent=True)
    except Exception:                                       # noqa: BLE001
        pass
    left, _, right = split_reaction(line) or ('', '', '')
    reactant, product = Chem.MolFromSmarts(left), Chem.MolFromSmarts(right)
    outcome: Dict[str, Any] = {
        'ok': True, 'level': 'ok',
        'reactants': built.GetNumReactantTemplates(),
        'products': built.GetNumProductTemplates(),
        'auto_maps': [], 'agents': '',
        'deleted_maps': [], 'deleted_unmapped': 0, 'new_atoms': 0,
    }
    if reactant is not None and product is not None:
        left_maps, right_maps = _maps(reactant), _maps(product)
        outcome['deleted_maps'] = [n for n in left_maps if n not in right_maps]
        outcome['deleted_unmapped'] = sum(1 for a in reactant.GetAtoms()
                                          if not a.GetAtomMapNum())
        outcome['new_atoms'] = sum(1 for a in product.GetAtoms()
                                   if not a.GetAtomMapNum())
        if not left_maps and not right_maps:
            outcome['level'] = 'note'
            outcome['status'] = ('Nothing is mapped -- the rule builds the '
                                 'product out of nothing.')
            return outcome
    outcome['status'] = describe(outcome)
    return outcome


def reaction_smarts_from_rxn_block(block: str) -> Optional[str]:
    """The fallback road, for when ``getSmarts`` will not write one.

    Ketcher's SMARTS export is not infallible -- a query feature on an
    aromatic ring has thrown rather than written -- and the version is
    deliberately not pinned, so there has to be a second way home.  A V3000
    RXN carries the atom maps and the aromatic bonds; plain element atoms come
    back as ``[#6]``, which against a ``:`` bond is the same query.
    """
    text = str(block or '')
    if not text.strip():
        return None
    try:
        drawn = rdChemReactions.ReactionFromRxnBlock(text, sanitize=False)
    except Exception:                                       # noqa: BLE001
        drawn = None
    if drawn is None:
        return None
    try:
        return rdChemReactions.ReactionToSmarts(drawn)
    except Exception:                                       # noqa: BLE001
        return None


def query_smarts_from_molblock(block: str) -> Optional[str]:
    """The same fallback for a pattern box, where there is no arrow."""
    text = str(block or '')
    if not text.strip():
        return None
    try:
        drawn = Chem.MolFromMolBlock(text, sanitize=False, removeHs=False)
    except Exception:                                       # noqa: BLE001
        drawn = None
    if drawn is None:
        return None
    try:
        return Chem.MolToSmarts(drawn)
    except Exception:                                       # noqa: BLE001
        return None


def rxn_block_from_smarts(line: str, *, v3000: bool = False) -> Optional[str]:
    """A rule as an RXN file, which is the only way its maps reach the editor.

    ``setMolecule`` reads a SMARTS happily and drops every atom map on the way
    in -- measured against Ketcher 3.17: ``[c&H1:1]:[c&H1:2]:...`` went in and
    ``[c;H]:[c;H]:c:c:c:c`` came back.  A rule opened that way and read again
    is a *different rule*, because the mapping then has to be guessed.  Sent
    as an RXN the maps survive, at the cost of the query features the atom
    block cannot state.
    """
    parts = split_reaction(line)
    if parts is None:
        return None
    left, _, right = parts
    reactant = Chem.MolFromSmarts(left)
    if reactant is None:
        return None
    # The product goes over as a real molecule, never as the query it is
    # written as.  A template bond that SMARTS leaves unspecified has no bond
    # order for the molfile writer to write, so it writes a single bond -- and
    # the editor then draws, and reads back, a saturated ring where an
    # aromatic one was meant.  ``concretize_product`` is what gives it real
    # aromatic bonds to write.
    drawn = concretize_product(right)
    product = Chem.MolFromSmiles(drawn, sanitize=False) if drawn else None
    if product is None:
        product = Chem.MolFromSmarts(right)
    if product is None:
        return None
    try:
        product.UpdatePropertyCache(strict=False)
    except Exception:                                       # noqa: BLE001
        pass
    try:
        built = rdChemReactions.ChemicalReaction()
        built.AddReactantTemplate(reactant)
        built.AddProductTemplate(product)
        return rdChemReactions.ReactionToRxnBlock(built, forceV3000=v3000)
    except Exception:                                       # noqa: BLE001
        return None


def survives_the_editor(line: str) -> bool:
    """Whether this rule comes back from a trip through Ketcher unchanged.

    Neither road into the editor is lossless.  A SMARTS keeps the query
    features and loses every atom map; an RXN keeps the maps and loses what
    the atom block cannot state -- Indigo ignores the SMARTS query sgroups
    RDKit writes into a V3000, so ``[c&H1:1]`` arrives as ``[#6:1]`` either
    way.  The maps are the half worth keeping, so the RXN is the road taken,
    and this says whether taking it costs anything for *this* rule, by making
    the trip and comparing what the two do rather than how they are spelled.

    A check against `_PROBES`, not a proof.  It answers "would this rule still
    do the same thing to these molecules", which is the question worth putting
    to somebody about to press the button; it can be wrong about a molecule
    that is not among them.
    """
    block = rxn_block_from_smarts(line)
    if not block:
        return False
    back = reaction_smarts_from_rxn_block(block)
    if not back:
        return False
    landed = normalize_reaction_smarts(back)
    if not landed.get('ok'):
        return False
    here = _behaviour(str(line), True)
    return here is not None and here == _behaviour(landed['smarts'], True)


def trial_on_seed(smarts: str, seed_smiles: str) -> Dict[str, Any]:
    """Run the rule against the seed that is actually in the box.

    "Reaction SMARTS produced no products" arrives after the Run button, names
    no rule and gives no reason, and the two things it usually means are very
    different: a template that never matched, and one that matched and built
    something RDKit will not accept. Both are cheap to find out at the moment
    the drawing is read, so they are said there.
    """
    seed = str(seed_smiles or '').strip()
    if not seed or '>>' not in str(smarts or ''):
        return {'level': '', 'status': ''}
    mol = Chem.MolFromSmiles(seed)
    if mol is None:
        mol = Chem.MolFromSmiles(seed, sanitize=False)
        if mol is not None:
            try:
                mol.UpdatePropertyCache(strict=False)
            except Exception:                               # noqa: BLE001
                mol = None
    if mol is None:
        return {'level': '', 'status': ''}
    try:
        rxn = rdChemReactions.ReactionFromSmarts(str(smarts))
    except Exception:                                       # noqa: BLE001
        rxn = None
    if rxn is None or not rxn.GetNumReactantTemplates():
        return {'level': '', 'status': ''}

    try:
        hits = len(mol.GetSubstructMatches(rxn.GetReactantTemplate(0)))
    except Exception:                                       # noqa: BLE001
        hits = 0
    if not hits:
        return {'level': 'note',
                'status': 'on the seed: no match — the bonds as drawn may be '
                          'single and double where the seed is aromatic'}

    made, refused = set(), ''
    try:
        groups = rxn.RunReactants((Chem.Mol(mol),))
    except Exception:                                       # noqa: BLE001
        groups = ()
    for group in groups:
        if not group:
            continue
        product = group[0]
        try:
            Chem.SanitizeMol(product)
            made.add(Chem.MolToSmiles(product))
        except Exception as exc:                            # noqa: BLE001
            refused = refused or str(exc).split('\n')[0][:90]
    if made:
        return {'level': 'ok',
                'status': f"on the seed: {len(made)} product"
                          f"{'s' if len(made) != 1 else ''}"}
    # The two ways a matching rule still builds nothing, in the order they
    # actually happen: an atom left unmapped is deleted, taking its ring bonds
    # with it, and an atom drawn aliphatic inside a ring the rest of which
    # stays aromatic cannot be sanitized.
    try:
        loose = sum(1 for atom in rxn.GetReactantTemplate(0).GetAtoms()
                    if not atom.GetAtomMapNum())
    except Exception:                                       # noqa: BLE001
        loose = 0
    if loose:
        advice = (f' — {loose} atom{"s" if loose != 1 else ""} on the left '
                  f'{"are" if loose != 1 else "is"} unmapped and so deleted; '
                  f'map {"them" if loose != 1 else "it"} if '
                  f'{"they" if loose != 1 else "it"} should stay')
    else:
        advice = ' — what stays aromatic has to be drawn aromatic'
    return {'level': 'note',
            'status': (f'on the seed: matches {hits}x, but no product survives'
                       + (f' — {refused}' if refused else '') + advice)}


def merge_maps(smarts: str, block: str) -> Tuple[str, int]:
    """Put back the atom maps Ketcher drops when it writes a SMARTS.

    Measured against Ketcher 3.17: a generic atom -- A, Q, X, M, the ones that
    say "anything here" -- loses its mapping in ``getSmarts`` and keeps it
    everywhere else.  ``A(:1)~C(:2)~A(:3)>>A(:1)~N(:2)~A(:3)`` comes back as
    ``[*;D3]~[#6:2;D2]~[*;D3]>>[*]~[#7:2]~[*]``, and a rule whose generic atoms
    are unmapped deletes and rebuilds them, so it matches and makes nothing.

    The same drawing's RXN carries all six maps, so they are taken from there.
    Only holes are filled, never an existing map overwritten, and only when
    the two readings agree on how many atoms each side has and on what they
    are -- otherwise the two orders are not the same order and nothing is
    touched.
    """
    # Repaired first: the string as Ketcher writes it does not parse, and an
    # atom carrying a query property is exactly where it puts the map wrong.
    text = repair_atom_maps(smarts)
    parts = split_reaction(text)
    if parts is None or not str(block or '').strip():
        return text, 0
    try:
        drawn = rdChemReactions.ReactionFromRxnBlock(str(block), sanitize=False)
    except Exception:                                       # noqa: BLE001
        drawn = None
    if drawn is None:
        return text, 0

    sides = [
        (parts[0], [drawn.GetReactantTemplate(i)
                    for i in range(drawn.GetNumReactantTemplates())]),
        (parts[2], [drawn.GetProductTemplate(i)
                    for i in range(drawn.GetNumProductTemplates())]),
    ]
    written, filled = [], 0
    for side, templates in sides:
        here = Chem.MolFromSmarts(side)
        if here is None or not templates:
            return text, 0
        # One mol per side either way: a '.' in the SMARTS is what several
        # templates are, so they line up when flattened in order.
        reference = templates[0]
        for extra in templates[1:]:
            reference = Chem.CombineMols(reference, extra)
        if here.GetNumAtoms() != reference.GetNumAtoms():
            return text, 0
        for mine, theirs in zip(here.GetAtoms(), reference.GetAtoms()):
            # Only where both readings name an element.  A custom query that
            # replaces the label -- "N or O" on a carbon -- is atomic number 0
            # in the SMARTS and 6 in the atom block, and it is still the same
            # atom in the same place; refusing that is refusing the case the
            # map is most often missing from.
            if not mine.GetAtomicNum() or not theirs.GetAtomicNum():
                continue
            if mine.GetAtomicNum() != theirs.GetAtomicNum():
                return text, 0
        for mine, theirs in zip(here.GetAtoms(), reference.GetAtoms()):
            if not mine.GetAtomMapNum() and theirs.GetAtomMapNum():
                mine.SetAtomMapNum(theirs.GetAtomMapNum())
                filled += 1
        try:
            written.append(Chem.MolToSmarts(here))
        except Exception:                                   # noqa: BLE001
            return text, 0
    if not filled:
        return text, 0
    return f'{written[0]}>>{written[1]}', filled
