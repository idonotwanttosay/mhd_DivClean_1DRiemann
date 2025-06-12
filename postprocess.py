import glob
import re
from typing import Tuple, List, Optional
import pandas as pd
import numpy as np


def available_steps(prefix: str) -> List[int]:
    """Return sorted step numbers for given prefix."""
    files = glob.glob(f"Result/out_{prefix}_*.csv")
    return sorted(int(re.findall(rf"_{prefix}_(\d+)", f)[0]) for f in files)


def read_grid(prefix: str = "rho") -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """Read grid coordinates from the first output file. Works for 1D or 2D."""
    f = glob.glob(f"Result/out_{prefix}_*.csv")[0]
    df = pd.read_csv(f, header=None)
    if df.shape[1] == 2:
        xs = df[0].values
        return xs, None
    xs = np.unique(df[0])
    ys = np.unique(df[1])
    return xs, ys


def load_field(prefix: str, step: int, xs: np.ndarray, ys: Optional[np.ndarray]):
    """Load a field and reshape according to the grid."""
    df = pd.read_csv(f"Result/out_{prefix}_{step}.csv", header=None)
    if ys is None:
        return df[1].values
    return df[2].values.reshape(len(xs), len(ys)).T

