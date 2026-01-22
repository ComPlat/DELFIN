"""Schema validation for CONTROL.txt configurations."""
from __future__ import annotations

from dataclasses import dataclass
import difflib
import re
from typing import Any, Callable, Iterable, Mapping, MutableMapping

ORCA_SOLVENTS = (
    "1,1,1-trichloroethane",
    "1,1,2-trichloroethane",
    "1,2,4-trimethylbenzene",
    "1,2-dibromoethane",
    "1,2-dichloroethane",
    "1,2-ethanediol",
    "1,4-dioxane",
    "dioxane",
    "1-bromo-2-methylpropane",
    "1-bromooctane",
    "bromooctane",
    "1-bromopentane",
    "1-bromopropane",
    "1-butanol",
    "butanol",
    "1-chlorohexane",
    "chlorohexane",
    "1-chloropentane",
    "1-chloropropane",
    "1-decanol",
    "decanol",
    "1-fluorooctane",
    "1-heptanol",
    "heptanol",
    "1-hexanol",
    "hexanol",
    "1-hexene",
    "1-hexyne",
    "1-iodobutane",
    "1-iodohexadecane",
    "hexadecyliodide",
    "1-iodopentane",
    "1-iodopropane",
    "1-nitropropane",
    "1-nonanol",
    "nonanol",
    "1-octanol",
    "octanol",
    "1-pentanol",
    "pentanol",
    "1-pentene",
    "1-propanol",
    "propanol",
    "2,2,2-trifluoroethanol",
    "2,2,4-trimethylpentane",
    "isooctane",
    "2,4-dimethylpentane",
    "2,4-dimethylpyridine",
    "2,6-dimethylpyridine",
    "2-bromopropane",
    "2-butanol",
    "secbutanol",
    "2-chlorobutane",
    "2-heptanone",
    "2-hexanone",
    "2-methoxyethanol",
    "methoxyethanol",
    "2-methyl-1-propanol",
    "isobutanol",
    "2-methyl-2-propanol",
    "2-methylpentane",
    "2-methylpyridine",
    "2methylpyridine",
    "2-nitropropane",
    "2-octanone",
    "2-pentanone",
    "2-propanol",
    "isopropanol",
    "2-propen-1-ol",
    "e-2-pentene",
    "3-methylpyridine",
    "3-pentanone",
    "4-heptanone",
    "4-methyl-2-pentanone",
    "4methyl2pentanone",
    "4-methylpyridine",
    "5-nonanone",
    "acetic acid",
    "aceticacid",
    "acetone",
    "acetonitrile",
    "mecn",
    "ch3cn",
    "acetophenone",
    "ammonia",
    "aniline",
    "anisole",
    "benzaldehyde",
    "benzene",
    "benzonitrile",
    "benzyl alcohol",
    "benzylalcohol",
    "bromobenzene",
    "bromoethane",
    "bromoform",
    "butanal",
    "butanoic acid",
    "butanone",
    "butanonitrile",
    "butyl ethanoate",
    "butyl acetate",
    "butylacetate",
    "butylamine",
    "n-butylbenzene",
    "butylbenzene",
    "sec-butylbenzene",
    "secbutylbenzene",
    "tert-butylbenzene",
    "tbutylbenzene",
    "carbon disulfide",
    "carbondisulfide",
    "cs2",
    "carbon tetrachloride",
    "ccl4",
    "chlorobenzene",
    "chloroform",
    "chcl3",
    "a-chlorotoluene",
    "o-chlorotoluene",
    "conductor",
    "m-cresol",
    "mcresol",
    "o-cresol",
    "cyclohexane",
    "cyclohexanone",
    "cyclopentane",
    "cyclopentanol",
    "cyclopentanone",
    "decalin",
    "cis-decalin",
    "n-decane",
    "decane",
    "dibromomethane",
    "dibutylether",
    "o-dichlorobenzene",
    "odichlorobenzene",
    "e-1,2-dichloroethene",
    "z-1,2-dichloroethene",
    "dichloromethane",
    "ch2cl2",
    "dcm",
    "diethyl ether",
    "diethylether",
    "diethyl sulfide",
    "diethylamine",
    "diiodomethane",
    "diisopropyl ether",
    "diisopropylether",
    "cis-1,2-dimethylcyclohexane",
    "dimethyl disulfide",
    "n,n-dimethylacetamide",
    "dimethylacetamide",
    "n,n-dimethylformamide",
    "dimethylformamide",
    "dmf",
    "dimethylsulfoxide",
    "dmso",
    "diphenylether",
    "dipropylamine",
    "n-dodecane",
    "dodecane",
    "ethanethiol",
    "ethanol",
    "ethyl acetate",
    "ethylacetate",
    "ethanoate",
    "ethyl methanoate",
    "ethyl phenyl ether",
    "ethoxybenzene",
    "ethylbenzene",
    "fluorobenzene",
    "formamide",
    "formic acid",
    "furan",
    "furane",
    "n-heptane",
    "heptane",
    "n-hexadecane",
    "hexadecane",
    "n-hexane",
    "hexane",
    "hexanoic acid",
    "iodobenzene",
    "iodoethane",
    "iodomethane",
    "isopropylbenzene",
    "p-isopropyltoluene",
    "isopropyltoluene",
    "mesitylene",
    "methanol",
    "methyl benzoate",
    "methyl butanoate",
    "methyl ethanoate",
    "methyl methanoate",
    "methyl propanoate",
    "n-methylaniline",
    "methylcyclohexane",
    "n-methylformamide",
    "methylformamide",
    "nitrobenzene",
    "phno2",
    "nitroethane",
    "nitromethane",
    "meno2",
    "o-nitrotoluene",
    "onitrotoluene",
    "n-nonane",
    "nonane",
    "n-octane",
    "octane",
    "n-pentadecane",
    "pentadecane",
    "octanol(wet)",
    "wetoctanol",
    "woctanol",
    "pentanal",
    "n-pentane",
    "pentane",
    "pentanoic acid",
    "pentyl ethanoate",
    "pentylamine",
    "perfluorobenzene",
    "hexafluorobenzene",
    "phenol",
    "propanal",
    "propanoic acid",
    "propanonitrile",
    "propyl ethanoate",
    "propylamine",
    "pyridine",
    "tetrachloroethene",
    "c2cl4",
    "tetrahydrofuran",
    "thf",
    "tetrahydrothiophene-s,s-dioxide",
    "tetrahydrothiophenedioxide",
    "sulfolane",
    "tetralin",
    "thiophene",
    "thiophenol",
    "toluene",
    "trans-decalin",
    "tributylphosphate",
    "trichloroethene",
    "triethylamine",
    "n-undecane",
    "undecane",
    "water",
    "h2o",
    "xylene",
    "m-xylene",
    "o-xylene",
    "p-xylene",
)

_SOLVENTS_LOWER = {name.lower(): name for name in ORCA_SOLVENTS}

ORCA_FUNCTIONALS = (
    "HFS",
    "LSD",
    "VWN5",
    "VWN3",
    "PWLDA",
    "BNULL",
    "BVWN",
    "BP",
    "BP86",
    "PW91",
    "mPWPW",
    "mPWLYP",
    "BLYP",
    "GP",
    "GLYP",
    "PBE",
    "revPBE",
    "RPBE",
    "PWP",
    "OLYP",
    "OPBE",
    "XLYP",
    "B97D",
    "PW86PBE",
    "RPW86PBE",
    "M06L",
    "TPSS",
    "revTPSS",
    "SCANfunc",
    "RSCAN",
    "R2SCAN",
    "B1LYP",
    "B1P",
    "G1LYP",
    "G1P",
    "B3LYP",
    "B3P",
    "G3LYP",
    "G3P",
    "PBE0",
    "PWP1",
    "mPW1PW",
    "mPW1LYP",
    "PW91_0",
    "O3LYP",
    "X3LYP",
    "B97",
    "BHANDHLYP",
    "TPSSh",
    "TPSS0",
    "PW6B95",
    "M06",
    "M062X",
    "r2SCANh",
    "r2SCAN0",
    "r2SCAN50",
    "wB97",
    "wB97X",
    "CAM-B3LYP",
    "LC_BLYP",
    "LC_PBE",
    "wr2SCAN",
)

_FUNCTIONALS_LOWER = {name.lower(): name for name in ORCA_FUNCTIONALS}

ORCA_BASIS_SETS = (
    "STO-3G",
    "MINI",
    "MINIS",
    "MINIX",
    "MIDI",
    "3-21G",
    "3-21GSP",
    "4-22GSP",
    "6-31G",
    "6-31G*",
    "m6-31G",
    "m6-31G*",
    "6-31G**",
    "6-31G(d)",
    "6-31G(d,p)",
    "6-31G(2d)",
    "6-31G(2d,p)",
    "6-31G(2d,2p)",
    "6-31G(2df)",
    "6-31G(2df,2p)",
    "6-31G(2df,2pd)",
    "6-31+G*",
    "6-31+G**",
    "6-31+G(d)",
    "6-31+G(d,p)",
    "6-31+G(2d)",
    "6-31+G(2d,p)",
    "6-31+G(2d,2p)",
    "6-31+G(2df)",
    "6-31+G(2df,2p)",
    "6-31+G(2df,2pd)",
    "6-31++G**",
    "6-31++G(d,p)",
    "6-31++G(2d,p)",
    "6-31++G(2d,2p)",
    "6-31++G(2df,2p)",
    "6-31++G(2df,2pd)",
    "6-311G",
    "6-311G*",
    "6-311G**",
    "6-311G(d)",
    "6-311G(d,p)",
    "6-311G(2d)",
    "6-311G(2d,p)",
    "6-311G(2d,2p)",
    "6-311G(2df)",
    "6-311G(2df,2p)",
    "6-311G(2df,2pd)",
    "6-311G(3df)",
    "6-311G(3df,3pd)",
    "6-311+G*",
    "6-311+G**",
    "6-311+G(d)",
    "6-311+G(d,p)",
    "6-311+G(2d)",
    "6-311+G(2d,p)",
    "6-311+G(2d,2p)",
    "6-311+G(2df)",
    "6-311+G(2df,2p)",
    "6-311+G(2df,2pd)",
    "6-311+G(3df)",
    "6-311+G(3df,2p)",
    "6-311+G(3df,3pd)",
    "6-311++G**",
    "6-311++G(d,p)",
    "6-311++G(2d,p)",
    "6-311++G(2d,2p)",
    "6-311++G(2df,2p)",
    "6-311++G(2df,2pd)",
    "6-311++G(3df,3pd)",
    "SV",
    "SV(P)",
    "SVP",
    "TZV",
    "TZV(P)",
    "TZVP",
    "TZVPP",
    "QZVP",
    "QZVPP",
    "DKH-SV(P)",
    "DKH-SVP",
    "DKH-TZV(P)",
    "DKH-TZVP",
    "DKH-TZVPP",
    "DKH-QZVP",
    "DKH-QZVPP",
    "ZORA-SV(P)",
    "ZORA-SVP",
    "ZORA-TZV(P)",
    "ZORA-TZVP",
    "ZORA-TZVPP",
    "ZORA-QZVP",
    "ZORA-QZVPP",
    "def2-mSVP",
    "def2-mTZVP",
    "def2-mTZVPP",
    "def2-SV(P)",
    "def2-SVP",
    "def2-TZVP(-f)",
    "def2-TZVP",
    "def2-TZVPP",
    "def2-QZVP",
    "def2-QZVPP",
    "def2-SVPD",
    "def2-TZVPD",
    "def2-TZVPPD",
    "def2-QZVPD",
    "def2-QZVPPD",
    "dhf-SV(P)",
    "dhf-SVP",
    "dhf-TZVP",
    "dhf-TZVPP",
    "dhf-QZVP",
    "dhf-QZVPP",
    "dhf-SV(P)-2c",
    "dhf-SVP-2c",
    "dhf-TZVP-2c",
    "dhf-TZVPP-2c",
    "dhf-QZVP-2c",
    "dhf-QZVPP-2c",
    "DKH-def2-SV(P)",
    "DKH-def2-SVP",
    "DKH-def2-TZVP(-f)",
    "DKH-def2-TZVP",
    "DKH-def2-TZVPP",
    "DKH-def2-QZVPP",
    "ZORA-def2-SV(P)",
    "ZORA-def2-SVP",
    "ZORA-def2-TZVP(-f)",
    "ZORA-def2-TZVP",
    "ZORA-def2-TZVPP",
    "ZORA-def2-QZVPP",
    "ma-def2-mSVP",
    "ma-def2-SV(P)",
    "ma-def2-SVP",
    "ma-def2-TZVP(-f)",
    "ma-def2-TZVP",
    "ma-def2-TZVPP",
    "ma-def2-QZVP",
    "ma-def2-QZVPP",
    "ma-DKH-def2-SV(P)",
    "ma-DKH-def2-SVP",
    "ma-DKH-def2-TZVP(-f)",
    "ma-DKH-def2-TZVP",
    "ma-DKH-def2-TZVPP",
    "ma-DKH-def2-QZVPP",
    "ma-ZORA-def2-SV(P)",
    "ma-ZORA-def2-SVP",
    "ma-ZORA-def2-TZVP(-f)",
    "ma-ZORA-def2-TZVP",
    "ma-ZORA-def2-TZVPP",
    "ma-ZORA-def2-QZVPP",
    "old-SV",
    "old-SV(P)",
    "old-SVP",
    "old-TZV",
    "old-TZV(P)",
    "old-TZVP",
    "old-TZVPP",
    "old-DKH-SV(P)",
    "old-DKH-SVP",
    "old-DKH-TZV(P)",
    "old-DKH-TZVP",
    "old-DKH-TZVPP",
    "old-ZORA-SV(P)",
    "old-ZORA-SVP",
    "old-ZORA-TZV(P)",
    "old-ZORA-TZVP",
    "old-ZORA-TZVPP",
    "ANO-SZ",
    "ANO-pVDZ",
    "ANO-pVTZ",
    "ANO-pVQZ",
    "ANO-pV5Z",
    "ANO-pV6Z",
    "aug-ANO-pVDZ",
    "aug-ANO-pVTZ",
    "aug-ANO-pVQZ",
    "aug-ANO-pV5Z",
    "saug-ANO-pVDZ",
    "saug-ANO-pVTZ",
    "saug-ANO-pVQZ",
    "saug-ANO-pV5Z",
    "ANO-RCC-DZP",
    "ANO-RCC-TZP",
    "ANO-RCC-QZP",
    "ANO-RCC-Full",
    "pc-0",
    "pc-1",
    "pc-2",
    "pc-3",
    "pc-4",
    "aug-pc-0",
    "aug-pc-1",
    "aug-pc-2",
    "aug-pc-3",
    "aug-pc-4",
    "pcJ-0",
    "pcJ-1",
    "pcJ-2",
    "pcJ-3",
    "pcJ-4",
    "aug-pcJ-0",
    "aug-pcJ-1",
    "aug-pcJ-2",
    "aug-pcJ-3",
    "aug-pcJ-4",
    "pcseg-0",
    "pcseg-1",
    "pcseg-2",
    "pcseg-3",
    "pcseg-4",
    "aug-pcseg-0",
    "aug-pcseg-1",
    "aug-pcseg-2",
    "aug-pcseg-3",
    "aug-pcseg-4",
    "pcSseg-0",
    "pcSseg-1",
    "pcSseg-2",
    "pcSseg-3",
    "pcSseg-4",
    "aug-pcSseg-0",
    "aug-pcSseg-1",
    "aug-pcSseg-2",
    "aug-pcSseg-3",
    "aug-pcSseg-4",
    "W1-mtsmall",
    "W1-DZ",
    "W1-TZ",
    "W1-QZ",
    "W1-Opt",
    "Sapporo-DZP-2012",
    "Sapporo-TZP-2012",
    "Sapporo-QZP-2012",
    "Sapporo-DKH3-DZP-2012",
    "Sapporo-DKH3-TZP-2012",
    "Sapporo-DKH3-QZP-2012",
    "LANL08",
    "LANL08(f)",
    "LANL2DZ",
    "LANL2TZ",
    "LANL2TZ(f)",
    "vDZP",
    "def-TZVP",
    "ma-def-TZVP",
    "HGBS-5",
    "HGBS-7",
    "HGBS-9",
    "HGBSP1-5",
    "HGBSP1-7",
    "HGBSP1-9",
    "HGBSP2-5",
    "HGBSP2-7",
    "HGBSP2-9",
    "HGBSP3-5",
    "HGBSP3-7",
    "HGBSP3-9",
    "AHGBS-5",
    "AHGBS-7",
    "AHGBS-9",
    "AHGBSP1-5",
    "AHGBSP1-7",
    "AHGBSP1-9",
    "AHGBSP2-5",
    "AHGBSP2-7",
    "AHGBSP2-9",
    "AHGBSP3-5",
    "AHGBSP3-7",
    "AHGBSP3-9",
    "cc-pVDZ",
    "cc-pVTZ",
    "cc-pVQZ",
    "cc-pV5Z",
    "cc-pV6Z",
    "aug-cc-pVDZ",
    "aug-cc-pVTZ",
    "aug-cc-pVQZ",
    "aug-cc-pV5Z",
    "aug-cc-pV6Z",
    "cc-pVD(+d)Z",
    "cc-pVT(+d)Z",
    "cc-pVQ(+d)Z",
    "cc-pV5(+d)Z",
    "apr-cc-pV(Q+d)Z",
    "may-cc-pV(T+d)Z",
    "may-cc-pV(Q+d)Z",
    "jun-cc-pV(D+d)Z",
    "jun-cc-pV(T+d)Z",
    "jun-cc-pV(Q+d)Z",
    "jul-cc-pV(D+d)Z",
    "jul-cc-pV(T+d)Z",
    "jul-cc-pV(Q+d)Z",
    "maug-cc-pV(D+d)Z",
    "maug-cc-pV(T+d)Z",
    "maug-cc-pV(Q+d)Z",
    "aug-cc-pVD(+d)Z",
    "aug-cc-pVT(+d)Z",
    "aug-cc-pVQ(+d)Z",
    "aug-cc-pV5(+d)Z",
    "aug-cc-pV6(+d)Z",
    "aug-cc-pVTZ-J",
    "cc-pCVDZ",
    "cc-pCVTZ",
    "cc-pCVQZ",
    "cc-pCV5Z",
    "cc-pCV6Z",
    "aug-cc-pCVDZ",
    "aug-cc-pCVTZ",
    "aug-cc-pCVQZ",
    "aug-cc-pCV5Z",
    "aug-cc-pCV6Z",
    "cc-pwCVDZ",
    "cc-pwCVTZ",
    "cc-pwCVQZ",
    "cc-pwCV5Z",
    "aug-cc-pwCVDZ",
    "aug-cc-pwCVTZ",
    "aug-cc-pwCVQZ",
    "aug-cc-pwCV5Z",
    "cc-pVDZ-PP",
    "cc-pVTZ-PP",
    "cc-pVQZ-PP",
    "cc-pV5Z-PP",
    "aug-cc-pVDZ-PP",
    "aug-cc-pVTZ-PP",
    "aug-cc-pVQZ-PP",
    "aug-cc-pV5Z-PP",
    "cc-pCVDZ-PP",
    "cc-pCVTZ-PP",
    "cc-pCVQZ-PP",
    "cc-pCV5Z-PP",
    "aug-cc-pCVDZ-PP",
    "aug-cc-pCVTZ-PP",
    "aug-cc-pCVQZ-PP",
    "aug-cc-pCV5Z-PP",
    "cc-pwCVDZ-PP",
    "cc-pwCVTZ-PP",
    "cc-pwCVQZ-PP",
    "cc-pwCV5Z-PP",
    "aug-cc-pwCVDZ-PP",
    "aug-cc-pwCVTZ-PP",
    "aug-cc-pwCVQZ-PP",
    "aug-cc-pwCV5Z-PP",
    "cc-pVDZ-DK",
    "cc-pVTZ-DK",
    "cc-pVQZ-DK",
    "cc-pV5Z-DK",
    "cc-pVDZ-DK3",
    "cc-pVTZ-DK3",
    "cc-pVQZ-DK3",
    "aug-cc-pVDZ-DK",
    "aug-cc-pVTZ-DK",
    "aug-cc-pVQZ-DK",
    "aug-cc-pV5Z-DK",
    "cc-pwCVDZ-DK",
    "cc-pwCVTZ-DK",
    "cc-pwCVQZ-DK",
    "cc-pwCV5Z-DK",
    "cc-pwCVDZ-DK3",
    "cc-pwCVTZ-DK3",
    "cc-pwCVQZ-DK3",
    "aug-cc-pwCVDZ-DK",
    "aug-cc-pwCVTZ-DK",
    "aug-cc-pwCVQZ-DK",
    "aug-cc-pwCV5Z-DK",
    "cc-pVDZ-F12",
    "cc-pVTZ-F12",
    "cc-pVQZ-F12",
    "cc-pVDZ-PP-F12",
    "cc-pVTZ-PP-F12",
    "cc-pVQZ-PP-F12",
    "cc-pCVDZ-F12",
    "cc-pCVTZ-F12",
    "cc-pCVQZ-F12",
    "haV(T+d)Z",
    "haV(Q+d)Z",
    "haV(5+d)Z",
    "Partridge-1",
    "Partridge-2",
    "Partridge-3",
    "Partridge-4",
    "x2c-SV(P)all",
    "x2c-SVPall",
    "x2c-TZVPall",
    "x2c-TZVPPall",
    "x2c-QZVPall",
    "x2c-QZVPPall",
    "x2c-SV(P)all-2c",
    "x2c-SVPall-2c",
    "x2c-TZVPall-2c",
    "x2c-TZVPPall-2c",
    "x2c-QZVPall-2c",
    "x2c-QZVPPall-2c",
    "x2c-SV(P)all-s",
    "x2c-SVPall-s",
    "x2c-TZVPall-s",
    "x2c-TZVPPall-s",
    "x2c-QZVPall-s",
    "x2c-QZVPPall-s",
    "x2c-QZVPall-2c-s",
    "x2c-QZVPPall-2c-s",
    "SARC-DKH-SVP",
    "SARC-DKH-TZVP",
    "SARC-DKH-TZVPP",
    "SARC-ZORA-SVP",
    "SARC-ZORA-TZVP",
    "SARC-ZORA-TZVPP",
    "SARC2-DKH-QZV",
    "SARC2-DKH-QZVP",
    "SARC2-ZORA-QZV",
    "SARC2-ZORA-QZVP",
    "D95",
    "D95p",
    "EPR-II",
    "EPR-III",
    "IGLO-II",
    "IGLO-III",
    "UGBS",
    "CP",
    "CP(PPP)",
    "Wachters+f",
    "def2/J",
    "def2-mTZVP/J",
    "def2-mTZVPP/J",
    "x2c/J",
    "SARC/J",
    "def2/JK",
    "def2/JKsmall",
    "cc-pVTZ/JK",
    "cc-pVQZ/JK",
    "cc-pV5Z/JK",
    "aug-cc-pVTZ/JK",
    "aug-cc-pVQZ/JK",
    "aug-cc-pV5Z/JK",
    "SARC2-DKH-QZV/JK",
    "SARC2-DKH-QZVP/JK",
    "SARC2-ZORA-QZV/JK",
    "SARC2-ZORA-QZVP/JK",
    "def2-SVP/C",
    "def2-TZVP/C",
    "def2-TZVPP/C",
    "def2-QZVPP/C",
    "def2-SVPD/C",
    "def2-TZVPD/C",
    "def2-TZVPPD/C",
    "def2-QZVPPD/C",
    "cc-pVDZ/C",
    "cc-pVTZ/C",
    "cc-pVQZ/C",
    "cc-pV5Z/C",
    "cc-pV6Z/C",
    "aug-cc-pVDZ/C",
    "aug-cc-pVTZ/C",
    "aug-cc-pVQZ/C",
    "aug-cc-pV5Z/C",
    "aug-cc-pV6Z/C",
    "cc-pwCVDZ/C",
    "cc-pwCVTZ/C",
    "cc-pwCVQZ/C",
    "cc-pwCV5Z/C",
    "aug-cc-pwCVDZ/C",
    "aug-cc-pwCVTZ/C",
    "aug-cc-pwCVQZ/C",
    "aug-cc-pwCV5Z/C",
    "cc-pVDZ-PP/C",
    "cc-pVTZ-PP/C",
    "cc-pVQZ-PP/C",
    "aug-cc-pVDZ-PP/C",
    "aug-cc-pVTZ-PP/C",
    "aug-cc-pVQZ-PP/C",
    "cc-pwCVDZ-PP/C",
    "cc-pwCVTZ-PP/C",
    "cc-pwCVQZ-PP/C",
    "aug-cc-pwCVDZ-PP/C",
    "aug-cc-pwCVTZ-PP/C",
    "aug-cc-pwCVQZ-PP/C",
    "cc-pVDZ-F12-MP2Fit",
    "cc-pVTZ-F12-MP2Fit",
    "cc-pVQZ-F12-MP2Fit",
    "cc-pVDZ-PP-F12-MP2Fit",
    "cc-pVTZ-PP-F12-MP2Fit",
    "cc-pVQZ-PP-F12-MP2Fit",
    "cc-pCVDZ-F12-MP2Fit",
    "cc-pCVTZ-F12-MP2Fit",
    "cc-pCVQZ-F12-MP2Fit",
    "cc-pVDZ-F12-CABS",
    "cc-pVTZ-F12-CABS",
    "cc-pVQZ-F12-CABS",
    "cc-pVDZ-F12-OptRI",
    "cc-pVTZ-F12-OptRI",
    "cc-pVQZ-F12-OptRI",
    "cc-pVDZ-PP-F12-OptRI",
    "cc-pVTZ-PP-F12-OptRI",
    "cc-pVQZ-PP-F12-OptRI",
    "aug-cc-pVDZ-PP-OptRI",
    "aug-cc-pVTZ-PP-OptRI",
    "aug-cc-pVQZ-PP-OptRI",
    "aug-cc-pV5Z-PP-OptRI",
    "cc-pCVDZ-F12-OptRI",
    "cc-pCVTZ-F12-OptRI",
    "cc-pCVQZ-F12-OptRI",
    "aug-cc-pwCVDZ-PP-OptRI",
    "aug-cc-pwCVTZ-PP-OptRI",
    "aug-cc-pwCVQZ-PP-OptRI",
    "aug-cc-pwCV5Z-PP-OptRI",
)

_BASIS_SETS_LOWER = {name.lower(): name for name in ORCA_BASIS_SETS}


def _suggest_from_options(value: str, options: Iterable[str]) -> list[str]:
    matches = difflib.get_close_matches(value, list(options), n=3, cutoff=0.6)
    return matches


@dataclass(frozen=True)
class FieldSpec:
    name: str
    coerce: Callable[[Any], Any]
    required: bool = False
    default: Any = None
    allow_none: bool = False


def _as_int(value: Any) -> int:
    if value is None or value == "":
        raise ValueError("must be an integer")
    return int(value)


def _as_float(value: Any) -> float:
    if value is None or value == "":
        raise ValueError("must be a float")
    return float(value)


def _as_non_negative_float(value: Any) -> float:
    if value is None or value == "":
        return 0.0
    parsed = float(value)
    if parsed < 0:
        raise ValueError("must be >= 0")
    return parsed


def _as_non_positive_float(value: Any) -> float:
    if value is None or value == "":
        return 0.0
    parsed = float(value)
    if parsed > 0:
        raise ValueError("must be <= 0 (imaginary frequencies are negative)")
    return parsed


def _as_str(value: Any) -> str:
    if value is None:
        return ""
    return str(value)


def _as_yes_no(value: Any) -> str:
    text = str(value or "no").strip().lower()
    return "yes" if text in {"yes", "true", "1", "on"} else "no"


def _as_charge(value: Any) -> int:
    if value is None or value == "":
        raise ValueError("must be an integer like -2, 0, or +3")
    text = str(value).strip()
    if not re.match(r"^[+-]?\d+$", text):
        raise ValueError("must be an integer like -2, 0, or +3")
    return int(text)


def _as_xtb_method(value: Any) -> str:
    text = str(value or "").strip().upper()
    if text in {"XTB0", "XTB1", "XTB2", "XTBFF"}:
        return text
    raise ValueError("must be XTB0, XTB1, XTB2, or XTBFF")


def _as_implicit_solvation_model(value: Any) -> str:
    text = str(value or "").strip().upper()
    if text == "":
        return ""
    if text == "C-PCM":
        return "CPCM"
    if text in {"CPCM", "SMD"}:
        return text
    raise ValueError("must be CPCM or SMD")


def _as_solvent(value: Any) -> str:
    text = str(value or "").strip()
    if text == "":
        return ""
    key = text.lower()
    if key in _SOLVENTS_LOWER:
        return _SOLVENTS_LOWER[key]
    suggestions = _suggest_from_options(key, _SOLVENTS_LOWER.keys())
    if suggestions:
        formatted = ", ".join(_SOLVENTS_LOWER[s] for s in suggestions)
        raise ValueError(f"unknown solvent '{text}'. Did you mean: {formatted}")
    raise ValueError(f"unknown solvent '{text}'")


def _as_functional(value: Any) -> str:
    text = str(value or "").strip()
    if text == "":
        raise ValueError("must be one of the ORCA functionals list")
    libxc_match = re.match(r"(?i)^libxc\((.+)\)$", text)
    if libxc_match:
        inner = libxc_match.group(1).strip()
        if not inner:
            raise ValueError("LibXC(...) requires a functional name inside parentheses")
        inner_key = inner.lower()
        if inner_key in _FUNCTIONALS_LOWER:
            canonical = _FUNCTIONALS_LOWER[inner_key]
            return f"LibXC({canonical})"
        return f"LibXC({inner})"
    key = text.lower()
    if key in _FUNCTIONALS_LOWER:
        return _FUNCTIONALS_LOWER[key]
    suggestions = _suggest_from_options(key, _FUNCTIONALS_LOWER.keys())
    if suggestions:
        formatted = ", ".join(_FUNCTIONALS_LOWER[s] for s in suggestions)
        raise ValueError(f"unknown functional '{text}'. Did you mean: {formatted}")
    raise ValueError(f"unknown functional '{text}'")


def _as_basis_set_optional(value: Any) -> str:
    text = str(value or "").strip()
    if text == "":
        return ""
    key = text.lower()
    if key in _BASIS_SETS_LOWER:
        return _BASIS_SETS_LOWER[key]
    suggestions = _suggest_from_options(key, _BASIS_SETS_LOWER.keys())
    if suggestions:
        formatted = ", ".join(_BASIS_SETS_LOWER[s] for s in suggestions)
        raise ValueError(f"unknown basis set '{text}'. Did you mean: {formatted}")
    raise ValueError(f"unknown basis set '{text}'")


def _as_basis_set_required(value: Any) -> str:
    text = str(value or "").strip()
    if text == "":
        raise ValueError("must be a valid ORCA basis set")
    return _as_basis_set_optional(text)


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, str):
        return [item.strip() for item in value.split(',') if item.strip()]
    raise ValueError("must be a list or comma-separated string")

def _as_parallel_strategy(value: Any) -> str:
    """Coerce user value into a known ORCA parallel strategy token."""
    text = str(value or "auto").strip().lower()
    if text in {"threads", "serial", "auto"}:
        return text
    if text in {"mpi", "default"}:
        return "auto"
    raise ValueError("must be one of: auto, threads, serial")

def _as_imag_scope(value: Any) -> str:
    """Coerce user value into a known IMAG scope."""
    text = str(value or "initial").strip().lower()
    if text in {"initial", "all"}:
        return text
    raise ValueError("must be one of: initial, all")

def _as_imag_option(value: Any) -> int:
    """Coerce IMAG scheduler behaviour selector."""
    if value is None or value == "":
        return 2
    try:
        option = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("must be 1 or 2") from exc
    if option not in (1, 2):
        raise ValueError("must be 1 or 2")
    return option


def _as_calc_potential_method(value: Any) -> int:
    """Coerce calc_potential_method: accepts 1, 2, or 3."""
    if value is None or value == "":
        return 2
    try:
        val = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("must be 1, 2, or 3") from exc
    if val not in (1, 2, 3):
        raise ValueError("must be 1, 2, or 3")
    return val


def _as_soscfhessup(value: Any) -> str:
    """Coerce deltaSCF_SOSCFHESSUP: accepts LSR1, LBFGS, LPOWELL, LBOFILL."""
    if value is None or value == "":
        return "LSR1"
    text = str(value).strip().upper()
    if text in {"LSR1", "LBFGS", "LPOWELL", "LBOFILL"}:
        return text
    raise ValueError("must be LSR1, LBFGS, LPOWELL, or LBOFILL")


_RELATIVITY_VALUES = {"ZORA", "X2C", "DKH", "DKH2"}


def _as_relativity(value: Any) -> str:
    """Coerce relativity: accepts ZORA, X2C, DKH (DKH2), or empty/none."""
    if value is None or value == "":
        return "none"
    text = str(value).strip()
    if not text:
        return "none"
    upper = text.upper()
    if upper in {"NONE", "NO"}:
        return "none"
    if upper in _RELATIVITY_VALUES:
        return upper
    raise ValueError("must be ZORA, X2C, DKH, DKH2, or none")


def _as_aux_jk(value: Any) -> str:
    """Coerce aux_jk: accepts def2/J, def2/JK, or empty."""
    if value is None or value == "":
        return ""
    text = str(value).strip()
    if not text:
        return ""
    upper = text.upper()
    if upper == "DEF2/J":
        return "def2/J"
    if upper == "DEF2/JK":
        return "def2/JK"
    raise ValueError("must be def2/J or def2/JK")


def _as_aux_jk_rel(value: Any) -> str:
    """Coerce aux_jk_rel: accepts SARC/J, x2c/J, cc-pVTZ/C, or empty."""
    if value is None or value == "":
        return ""
    text = str(value).strip()
    if not text:
        return ""
    upper = text.upper()
    if upper == "SARC/J":
        return "SARC/J"
    if upper == "X2C/J":
        return "x2c/J"
    if upper == "CC-PVTZ/C":
        return "cc-pVTZ/C"
    raise ValueError("must be SARC/J, x2c/J, or cc-pVTZ/C")


_GEOM_OPT_BASES = {
    "OPT",
    "TIGHTOPT",
    "VERYTIGHTOPT",
}


def _as_geom_opt(value: Any) -> str:
    """Coerce geom_opt tokens to allowed ORCA geometry optimization keywords."""
    if value is None or value == "":
        return "OPT"
    text = str(value).strip()
    if not text:
        return "OPT"
    normalized = " ".join(text.split()).upper()
    tokens = normalized.split()
    if tokens and tokens[0] in _GEOM_OPT_BASES:
        # Allow additional ORCA keywords like NOTRAH, NODIIS, etc.
        return normalized
    raise ValueError(
        "must start with OPT, TIGHTOPT, or VERYTIGHTOPT (additional keywords allowed)"
    )


def _as_freq_type(value: Any) -> str:
    """Coerce freq_type: accepts FREQ or numFREQ."""
    if value is None or value == "":
        return "FREQ"
    text = str(value).strip()
    if not text:
        return "FREQ"
    upper = text.upper()
    if upper in {"FREQ", "NUMFREQ"}:
        return "numFREQ" if upper == "NUMFREQ" else "FREQ"
    raise ValueError("must be FREQ or numFREQ")


def _as_initial_guess(value: Any) -> str:
    """Coerce initial_guess: accepts PAtom, PModel, Hueckel, or HCore."""
    if value is None or value == "":
        return "PModel"
    text = str(value).strip()
    if not text:
        return "PModel"
    upper = text.upper()
    if upper in {"PATOM", "PMODEL", "HUECKEL", "HCORE"}:
        if upper == "PATOM":
            return "PAtom"
        if upper == "PMODEL":
            return "PModel"
        if upper == "HUECKEL":
            return "Hueckel"
        return "HCore"
    raise ValueError("must be PAtom, PModel, Hueckel, or HCore")


_RI_JKX_KEYWORDS = {
    "RI",
    "NORI",
    "SPLITRIJ",
    "NOSPLITRIJ",
    "NOSPLITRIJK",
    "RIJONX",
    "RIJDX",
    "RIJCOSX",
    "NORIJCOSX",
    "NOCOSX",
    "RIJK",
    "NOSFITTING",
}


def _as_ri_jkx(value: Any) -> str:
    """Coerce ri_jkx: accepts ORCA RI/SplitRIJ/RIJCOSX/RIJK style keywords or empty."""
    if value is None or value == "":
        return ""
    text = str(value).strip()
    normalized = re.sub(r"[-_\\s]", "", text).upper()
    if normalized in _RI_JKX_KEYWORDS:
        return text
    raise ValueError(
        "must be one of: RI, NORI, SplitRIJ, NoSplitRIJ/NoSplitRIJK, RIJONX, RIJDX, "
        "RIJCOSX, NORIJCOSX, NOCOSX, RIJK/RI-JK, NoSFitting"
    )


# Dispersion correction: valid values and functional compatibility
DISP_CORR_VALUES = {"D2", "D3", "D3BJ", "D3ZERO", "D30", "D3TZ", "D4", "ABC", "ATM", "NOVDW", ""}

# Functionals compatible with D4
D4_FUNCTIONALS = {
    "HF", "BLYP", "BPBE", "BP86", "BPW91", "GLYP", "LB94", "MPWLYP", "MPWPW", "OLYP",
    "OPBE", "PBE", "RPBE", "REVPBE", "PW86PBE", "RPW86PBE", "PW91", "PW91P86", "XLYP",
    "B97BECKE", "TPSS", "REVTPSS", "SCAN", "B1LYP", "B3LYP", "BHANDHLYP", "B1P", "B3P86",
    "B1PW91", "B3PW91", "O3LYP", "REVPBE0", "REVPBE38", "PBE0", "PWP1", "PW1PW", "MPW1PW",
    "MPW1LYP", "PW6B95", "TPSSH", "TPSS0", "X3LYP", "M06L", "M06", "WB97", "WB97X",
    "B97M-D4", "CAM-B3LYP", "LC-BLYP", "B2PLYP", "B2GP-PLYP", "MPW2PLYP", "PWPB95",
    "B97-D", "RSCAN", "R2SCAN", "R2SCANH", "R2SCAN0", "R2SCAN50", "WB97X-D4REV",
    "WB97M-D4REV", "WR2SCAN", "R2SCAN0-DH", "R2SCAN-CIDH", "R2SCAN-QIDH", "R2SCAN0-2",
    "PR2SCAN50", "KPR2SCAN50", "WPR2SCAN50", "PR2SCAN69", "REVDSD-PBEP86/2021",
    "REVDOD-PBEP86/2021", "LRC-PBE",
}

# Functionals compatible with D3BJ
D3BJ_FUNCTIONALS = {
    "HF", "BP86", "BLYP", "REVPBE", "B97-D", "PBE", "RPBE", "RPW86PBE", "B3LYP",
    "BHANDHLYP", "TPSS", "TPSS0", "PBE0", "REVPBE38", "PW6B95", "B2PLYP", "MPWLYP",
    "OLYP", "BPBE", "OPBE", "B3PW91", "REVPBE0", "TPSSH", "CAM-B3LYP", "B2GP-PLYP",
    "PWPB95", "SCAN", "RSCAN", "R2SCAN", "R2SCANH", "R2SCAN0", "R2SCAN50", "WR2SCAN",
    "R2SCAN0-DH", "R2SCAN-CIDH", "R2SCAN-QIDH", "R2SCAN0-2", "PR2SCAN50", "KPR2SCAN50",
    "WPR2SCAN50", "PR2SCAN69", "REVDSD-PBEP86/2021", "REVDOD-PBEP86/2021", "DSD-BLYP",
    "DSD-BLYP/2013", "DSD-PBEB95", "DSD-PBEP86", "DSD-PBEP86/2013", "B97M-D3BJ",
    "WB97X-D3BJ", "WB97M-D3BJ", "WB97X-2", "PBE0DH", "PBE02", "PBE-QIDH",
}

# Functionals compatible with D3ZERO (D3(0), D3)
D3ZERO_FUNCTIONALS = {
    "HF", "BLYP", "BP86", "B97-D", "REVPBE", "PBE", "RPBE", "TPSS", "B3LYP", "PBE0",
    "PW6B95", "TPSS0", "B2PLYP", "B2GP-PLYP", "PWPB95", "MPWLYP", "BPBE", "BHANDHLYP",
    "TPSSH", "REVPBE0", "REVPBE38", "RPW86PBE", "B3PW91", "M06L", "M06", "M062X",
    "WB97X-D3", "CAM-B3LYP", "SCAN", "WB97X-2", "PBE0DH", "PBE02", "PBE-QIDH",
}


def _as_disp_corr(value: Any) -> str:
    """Coerce disp_corr: accepts D2, D3, D3BJ, D3ZERO, D30, D3TZ, D4, ABC, ATM, NOVDW or empty."""
    if value is None or value == "":
        return ""
    text = str(value).strip().upper()
    if text in DISP_CORR_VALUES:
        return text
    raise ValueError("must be D2, D3, D3BJ, D3ZERO, D30, D3TZ, D4, ABC, ATM, NOVDW, or empty")


def validate_disp_corr_functional_combo(disp_corr: str, functional: str) -> list[str]:
    """Validate dispersion correction and functional combination. Returns list of warnings."""
    warnings = []
    if not disp_corr:
        return warnings

    func_text = functional.strip()
    libxc_match = re.match(r"(?i)^libxc\((.+)\)$", func_text)
    if libxc_match:
        func_text = libxc_match.group(1).strip()
    func_upper = func_text.upper().replace("_", "-")
    disp_upper = disp_corr.upper()

    if disp_upper == "D4":
        if func_upper not in D4_FUNCTIONALS:
            warnings.append(f"D4 dispersion may not be parametrized for functional '{functional}'")
    elif disp_upper in {"D3BJ", "D3"}:
        if func_upper not in D3BJ_FUNCTIONALS:
            warnings.append(f"D3BJ dispersion may not be parametrized for functional '{functional}'")
    elif disp_upper in {"D3ZERO", "D30"}:
        if func_upper not in D3ZERO_FUNCTIONALS:
            warnings.append(f"D3ZERO dispersion may not be parametrized for functional '{functional}'")

    return warnings


def _as_method(value: Any) -> str:
    text = str(value or "").strip().lower()
    if text in {"classic", "manually", "occupier"}:
        return text
    raise ValueError("must be one of: classic, manually, OCCUPIER")


def _as_esd_modus(value: Any) -> str:
    text = str(value or "").strip().lower()
    if text in {"tddft", "deltascf", "hybrid1"}:
        return text
    raise ValueError("must be one of: TDDFT, deltaSCF, hybrid1")


def _as_properties_of_interest(value: Any) -> list[str] | str:
    if value is None or value == "":
        return ""
    if isinstance(value, (list, tuple)):
        items = [str(item).strip() for item in value if str(item).strip()]
    else:
        text = str(value).strip()
        if not text:
            return ""
        text = text.strip("[]").replace("'", "").replace('"', '')
        items = [item.strip() for item in text.split(",") if item.strip()]
    normalized = [item.upper() for item in items]
    for item in normalized:
        if item not in {"IP", "EA"}:
            raise ValueError("must be IP, EA, both, or empty")
    seen = set()
    unique = []
    for item in normalized:
        if item in seen:
            continue
        seen.add(item)
        unique.append(item)
    return unique


def _as_ics(value: Any) -> list[str] | str:
    if value is None or value == "":
        return ""
    if isinstance(value, (list, tuple)):
        items = [str(item).strip() for item in value if str(item).strip()]
    else:
        text = str(value).strip()
        if not text:
            return ""
        text = text.strip("[]").replace("'", "").replace('"', '')
        items = [item.strip() for item in text.split(",") if item.strip()]
    normalized = []
    for item in items:
        if not re.match(r"^[ST]\d+>(S1|T1)$", item):
            raise ValueError("ICs must be transitions like S2>S1 or T2>T1")
        normalized.append(item)
    return normalized


def _as_states(value: Any) -> list[str] | str:
    if value is None or value == "":
        return ""
    if isinstance(value, (list, tuple)):
        items = [str(item).strip() for item in value if str(item).strip()]
    else:
        text = str(value).strip()
        if not text:
            return ""
        text = text.strip("[]").replace("'", "").replace('"', '')
        items = [item.strip() for item in text.split(",") if item.strip()]
    normalized = []
    for item in items:
        match = re.match(r"^([ST])(\d+)$", item)
        if not match:
            raise ValueError("states must be S1..S6 or T1..T6")
        idx = int(match.group(2))
        if idx < 1 or idx > 6:
            raise ValueError("states must be S1..S6 or T1..T6")
        normalized.append(f"{match.group(1)}{idx}")
    return normalized


def _as_occupier_method(value: Any) -> str:
    text = str(value or "auto").strip().lower()
    if text in {"manual", "manually"}:
        return "manually"
    if text == "auto":
        return "auto"
    raise ValueError("must be auto or manually")


def _as_occupier_tree(value: Any) -> str:
    """Coerce user value into a known tree token (flat, deep2, deep3, deep, own)."""
    # If empty/None, use default "deep"
    if not value or str(value).strip() == "":
        return "deep"

    text = str(value).strip().lower()
    if text in {"deep", "tree"}:
        return "deep"
    if text in {"flat", "flatt", "legacy"}:
        return "flat"
    if text == "deep2":
        return "deep2"
    if text == "deep3":
        return "deep3"
    if text in {"deep4", "dee4"}:
        return "own"
    if text in {"own", "custom"}:
        return "own"
    # Legacy aliases for backwards compatibility - map to "deep"
    if text in {"deep5", "dee5", "deep6", "dee6"}:
        return "deep"
    # Default to "own" for invalid/unrecognized values (including multi-value syntax like "deep|flat|own")
    return "own"


def _as_ap_method(value: Any) -> int | None:
    """Coerce APMethod tokens into ORCA-compatible integers."""
    if value is None:
        return None

    text = str(value).strip()
    lowered = text.lower()

    if not text:
        return None
    if lowered in {"none", "off", "disable", "disabled", "false"}:
        return None
    if lowered == "apbs":
        return 2
    candidate = text
    if lowered.startswith("ap") and lowered[2:].isdigit():
        candidate = lowered[2:]

    try:
        parsed = int(candidate)
    except Exception as exc:  # noqa: BLE001
        raise ValueError("must be an integer (e.g. 1-4) or 'none'") from exc

    if parsed <= 0:
        return None
    if parsed not in {1, 2, 3}:
        raise ValueError("must be 1, 2, 3, or none")
    return parsed


CONTROL_FIELD_SPECS: Iterable[FieldSpec] = (
    FieldSpec("NAME", _as_str, default=""),
    FieldSpec("SMILES", _as_str, default=""),
    FieldSpec("charge", _as_charge, required=True),
    FieldSpec("multiplicity_global_opt", _as_int, allow_none=True),
    FieldSpec("PAL", _as_int, default=6),
    FieldSpec("number_explicit_solv_molecules", _as_int, default=0),
    FieldSpec("method", _as_method, required=True),
    FieldSpec("frequency_calculation", _as_yes_no, default="no"),
    FieldSpec("frequency_calculation_OCCUPIER", _as_yes_no, default="no"),
    FieldSpec("xTB_method", _as_xtb_method, default="XTB2"),
    FieldSpec("implicit_solvation_model", _as_implicit_solvation_model, default="CPCM"),
    FieldSpec("solvent", _as_solvent, default=""),
    FieldSpec("functional", _as_functional, default="PBE0"),
    FieldSpec("ri_jkx", _as_ri_jkx, default="RIJCOSX"),
    FieldSpec("aux_jk", _as_aux_jk, default="def2/J"),
    FieldSpec("aux_jk_rel", _as_aux_jk_rel, default="SARC/J"),
    FieldSpec("main_basisset", _as_basis_set_required, default="def2-SVP"),
    FieldSpec("main_basisset_rel", _as_basis_set_optional, default=""),
    FieldSpec("metal_basisset", _as_basis_set_optional, default=""),
    FieldSpec("metal_basisset_rel", _as_basis_set_optional, default=""),
    FieldSpec("initial_guess", _as_initial_guess, default="PModel"),
    FieldSpec("relativity", _as_relativity, default="none"),
    FieldSpec("geom_opt", _as_geom_opt, default="OPT"),
    FieldSpec("geom_opt_OCCUPIER", _as_geom_opt, default="OPT"),
    FieldSpec("freq_type", _as_freq_type, default="FREQ"),
    FieldSpec("orca_parallel_strategy", _as_parallel_strategy, default="auto"),
    FieldSpec("IMAG_scope", _as_imag_scope, default="initial"),
    FieldSpec("IMAG_option", _as_imag_option, default=2),
    FieldSpec("allow_imaginary_freq", _as_non_positive_float, default=0.0),
    FieldSpec("calc_potential_method", _as_calc_potential_method, default=2),
    FieldSpec("deltaSCF_SOSCFHESSUP", _as_soscfhessup, default="LSR1"),
    FieldSpec("OCCUPIER_method", _as_occupier_method, default="auto"),
    FieldSpec("OCCUPIER_tree", _as_occupier_tree, default="deep"),
    FieldSpec("OWN_progressive_from", _as_yes_no, default="no"),
    FieldSpec("OWN_TREE_PURE_WINDOW", _as_int, default=None, allow_none=True),
    FieldSpec("approximate_spin_projection_APMethod", _as_ap_method, default=2),
    FieldSpec("ESD_modus", _as_esd_modus, default="tddft"),
    FieldSpec("ESD_nroots", _as_int, default=15),
    FieldSpec("ESD_maxdim", _as_int, default=None, allow_none=True),
    FieldSpec("ESD_SOC", _as_yes_no, default="false"),
    FieldSpec("properties_of_interest", _as_properties_of_interest, default=""),
    FieldSpec("ICs", _as_ics, default=""),
    FieldSpec("states", _as_states, default=""),
)


def validate_control_config(config: MutableMapping[str, Any]) -> dict[str, Any]:
    """Validate and coerce CONTROL configuration values."""
    errors: list[str] = []
    validated: dict[str, Any] = dict(config)

    if "OCCUPIER_method" in config and str(config.get("OCCUPIER_method", "")).strip() == "":
        errors.append("OCCUPIER_method must be auto or manually")

    # First pass: validate OCCUPIER_method to determine if tree validation is needed
    occupier_method_raw = config.get("OCCUPIER_method", None)
    occupier_method = "manually"  # default
    if occupier_method_raw is not None and occupier_method_raw != "":
        try:
            occupier_method = _as_occupier_method(occupier_method_raw)
        except Exception:  # noqa: BLE001
            occupier_method = "manually"

    # Check if ESD_modul is enabled
    esd_modul_raw = config.get("ESD_modul", "no")
    esd_modul_enabled = str(esd_modul_raw).strip().lower() == "yes"

    for spec in CONTROL_FIELD_SPECS:
        # Skip OCCUPIER_tree validation if OCCUPIER_method is 'manually'
        if spec.name == "OCCUPIER_tree" and occupier_method == "manually":
            # Just set default without validation when manually mode
            validated[spec.name] = "deep"
            continue

        # Skip ESD_modus validation if ESD_modul is not enabled
        if spec.name == "ESD_modus" and not esd_modul_enabled:
            validated[spec.name] = spec.default
            continue

        raw = config.get(spec.name, None)
        if raw is None or raw == "":
            if spec.required and spec.default is None:
                errors.append(f"Missing required key: {spec.name}")
                continue
            if spec.default is not None:
                validated[spec.name] = spec.default
                continue
            if spec.allow_none:
                validated[spec.name] = None
                continue
        try:
            validated[spec.name] = spec.coerce(raw)
        except Exception as exc:  # noqa: BLE001
            errors.append(f"Invalid value for {spec.name}: {exc}")

    def _validate_sequence_list(seq_value: Any, label: str) -> None:
        if not isinstance(seq_value, list):
            errors.append(f"{label} must be a list of mappings")
            return
        for idx, item in enumerate(seq_value, start=1):
            if not isinstance(item, Mapping):
                errors.append(f"{label}[{idx}] must be a mapping")
                continue
            if "index" not in item or "m" not in item:
                errors.append(f"{label}[{idx}] must define 'index' and 'm'")
                continue
            try:
                int(item["index"])
                int(item["m"])
            except Exception:  # noqa: BLE001
                errors.append(f"{label}[{idx}] has non-integer 'index' or 'm'")

    # ensure electron sequences have expected structure if present
    for seq_key in ("even_seq", "odd_seq"):
        if seq_key in config:
            _validate_sequence_list(config[seq_key], seq_key)

    blocks = config.get("_occupier_sequence_blocks")
    if blocks is not None:
        if not isinstance(blocks, list):
            errors.append("_occupier_sequence_blocks must be a list")
        else:
            for block_idx, block in enumerate(blocks, start=1):
                if not isinstance(block, Mapping):
                    errors.append(f"sequence block #{block_idx} must be a mapping")
                    continue
                deltas = block.get("deltas")
                if not isinstance(deltas, list) or not all(isinstance(d, int) for d in deltas):
                    errors.append(f"sequence block #{block_idx} has invalid 'deltas'")
                if "even_seq" in block:
                    _validate_sequence_list(block["even_seq"], f"sequence block #{block_idx} even_seq")
                if "odd_seq" in block:
                    _validate_sequence_list(block["odd_seq"], f"sequence block #{block_idx} odd_seq")

    # Also report dispersion/functional compatibility if both fields parse.
    disp_raw = config.get("disp_corr", "")
    func_raw = config.get("functional", "")
    disp_val = None
    func_val = None
    try:
        disp_val = _as_disp_corr(disp_raw)
    except Exception as exc:  # noqa: BLE001
        errors.append(f"disp_corr invalid: {disp_raw!r} ({exc})")
    try:
        func_val = _as_functional(func_raw)
    except Exception:  # noqa: BLE001
        func_val = None
    if disp_val is not None and func_val is not None:
        errors.extend(validate_disp_corr_functional_combo(disp_val, func_val))

    if errors:
        raise ValueError("; ".join(errors))

    return validated
