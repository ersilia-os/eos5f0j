# imports
import os
import csv
import sys
import shutil
import platform
import subprocess
import tempfile
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from ersilia_pack_utils.core import read_smiles, write_out

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# canonical functional group table (committed next to this script)
FG_FILE = os.path.join(root, "functional_groups.csv")


def find_checkmol():
    """Locate the checkmol binary.

    Prebuilt, FP-exception-masked binaries for the two Linux architectures Ersilia
    targets are committed next to this script (checkmol-linux-x86_64 /
    checkmol-linux-aarch64). We resolve the architecture-specific binary directly so
    the model works in the packed Docker image regardless of whether
    install_checkmol.sh ran during the build. Fallbacks: a plain `checkmol` next to
    this script (e.g. compiled on a dev machine), the interpreter's bin directory,
    and PATH.
    """
    candidates = []
    if sys.platform.startswith("linux"):
        arch = {
            "x86_64": "x86_64",
            "amd64": "x86_64",
            "aarch64": "aarch64",
            "arm64": "aarch64",
        }.get(platform.machine().lower())
        if arch:
            candidates.append(os.path.join(root, "checkmol-linux-" + arch))
    candidates += [
        os.path.join(root, "checkmol"),
        os.path.join(os.path.dirname(sys.executable), "checkmol"),
    ]
    for path in candidates:
        if os.path.isfile(path):
            if not os.access(path, os.X_OK):
                # the executable bit can be lost in transit; restore it best-effort
                try:
                    os.chmod(path, 0o755)
                except OSError:
                    pass
            if os.access(path, os.X_OK):
                return path
    on_path = shutil.which("checkmol")
    if on_path:
        return on_path
    raise FileNotFoundError(
        "checkmol binary not found. Prebuilt binaries are committed next to this "
        "script (checkmol-linux-*); install_checkmol.sh installs/compiles them."
    )


CHECKMOL = find_checkmol()


def load_functional_groups():
    """Load the fixed, ordered list of checkmol functional groups.

    Returns the ordered list of feature names (one per functional group) and a
    mapping from checkmol's 1-based functional group number to the column index.
    """
    features = []
    fgnum_to_idx = {}
    with open(FG_FILE) as f:
        reader = csv.DictReader(f)
        for idx, row in enumerate(reader):
            features.append(row["feature"])
            fgnum_to_idx[int(row["fg_number"])] = idx
    return features, fgnum_to_idx


def smiles_to_molblock(smiles):
    """Convert a SMILES string to an MDL molblock; returns None if unparseable."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    # checkmol consumes MDL molfiles; 2D coordinates avoid degenerate geometry
    AllChem.Compute2DCoords(mol)
    return Chem.MolToMolBlock(mol)


def run_checkmol(molblock):
    """Run checkmol in position mode (-p) on a molblock and return {fg_number: count}.

    checkmol -p prints one line per detected functional group in the form
    '#NNN:count:positions' (no label), where NNN is the 1-based functional group
    number and count is the number of occurrences.
    """
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".mol", delete=False
        ) as tmp:
            tmp.write(molblock)
            tmp_path = tmp.name
        proc = subprocess.run(
            [CHECKMOL, "-p", tmp_path],
            capture_output=True,
            text=True,
            timeout=120,
        )
    finally:
        if tmp_path is not None and os.path.exists(tmp_path):
            os.remove(tmp_path)

    counts = {}
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line.startswith("#"):
            continue
        parts = line.split(":")
        if len(parts) < 2:
            continue
        try:
            fgnum = int(parts[0][1:])
            count = int(parts[1])
        except ValueError:
            continue
        counts[fgnum] = count
    return counts


def checkmol_wrapper(smiles_list, features, fgnum_to_idx):
    n_feat = len(features)
    results = []
    for smiles in smiles_list:
        molblock = smiles_to_molblock(smiles)
        if molblock is None:
            # unparseable SMILES -> row of None (empty cell in csv, NaN in bin)
            results.append([None] * n_feat)
            continue
        counts = run_checkmol(molblock)
        # integer counts (Python ints) so csv output is written as 0/1, not 0.0/1.0
        vector = [0] * n_feat
        for fgnum, count in counts.items():
            idx = fgnum_to_idx.get(fgnum)
            if idx is not None:
                vector[idx] = count
        results.append(vector)
    return results


# read SMILES from .csv (or .bin) file
_, smiles_list = read_smiles(input_file)

# load the fixed functional group schema
features, fgnum_to_idx = load_functional_groups()

# run model
outputs = checkmol_wrapper(smiles_list, features, fgnum_to_idx)

# check input and output have the same length
assert len(smiles_list) == len(outputs)

# write output. Values are integer counts; csv keeps them as integers and the
# binary path uses float32 so any failed compound can carry NaN.
write_out(outputs, features, output_file, np.float32)
