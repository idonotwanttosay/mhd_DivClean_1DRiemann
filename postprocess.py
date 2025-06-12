import glob
import os
import re
from typing import Tuple, List, Optional
import pandas as pd
import numpy as np


def available_steps(prefix: str) -> List[int]:
    """Return sorted step numbers for given prefix."""
    files = glob.glob(f"Result/out_{prefix}_*.csv")
    pattern = rf"_{prefix}_(\d+)"
    if not files:
        # Fallback to consolidated state files
        files = glob.glob("Result/out_state_*.csv")
        pattern = r"_state_(\d+)"
    return sorted(int(re.findall(pattern, f)[0]) for f in files)


def read_grid(prefix: str = "rho") -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """Read grid coordinates from the first output file. Works for 1D or 2D."""
    files = glob.glob(f"Result/out_{prefix}_*.csv")
    if files:
        f = files[0]
        df = pd.read_csv(f, header=None)
        if df.shape[1] == 2:
            xs = df[0].values
            return xs, None
        xs = np.unique(df[0])
        ys = np.unique(df[1])
        return xs, ys
    # Fallback to state files (assume 1D)
    f = glob.glob("Result/out_state_*.csv")[0]
    df = pd.read_csv(f)
    xs = df['x'].values
    return xs, None


def load_field(prefix: str, step: int, xs: np.ndarray, ys: Optional[np.ndarray]):
    """Load a field and reshape according to the grid."""
    fname = f"Result/out_{prefix}_{step}.csv"
    if os.path.exists(fname):
        df = pd.read_csv(fname, header=None)
        if ys is None:
            return df[1].values
        return df[2].values.reshape(len(xs), len(ys)).T
    # Fall back to consolidated state file
    fname = f"Result/out_state_{step}.csv"
    df = pd.read_csv(fname)
    return df[prefix].values

