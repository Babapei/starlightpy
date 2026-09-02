"""Non-negative regularization of x_j. Never Metropolis / annealing."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import nnls

# Couple neighboring ages; 3 dex young-vs-old is exp(-6) ~ 0 and is ignored.
_AGE_DEX_SCALE = 0.5
_MIN_COUPLE = 0.05


def _validate_ages(ages: ArrayLike, n_comp: int) -> NDArray[np.float64]:
    a = np.asarray(ages, dtype=float)
    if a.shape != (n_comp,):
        raise ValueError("template_ages must have one value per template.")
    if np.any(a <= 0) or not np.all(np.isfinite(a)):
        raise ValueError("template_ages must be finite and positive.")
    return a


def age_bin_ids(ages: NDArray[np.float64], edges: Optional[Sequence[float]]) -> NDArray[np.int_]:
    if edges is None:
        uniq = np.unique(ages)
        ids = np.empty(ages.size, dtype=int)
        for i, u in enumerate(uniq):
            ids[np.isclose(ages, u)] = i
        return ids
    edge_arr = np.asarray(edges, dtype=float)
    if edge_arr.size < 2:
        raise ValueError("age_bin_edges must contain at least two edges.")
    if np.any(np.diff(edge_arr) <= 0):
        raise ValueError("age_bin_edges must be strictly increasing.")
    ids = np.digitize(ages, edge_arr[1:-1], right=False)
    below = ages < edge_arr[0]
    above = ages >= edge_arr[-1]
    if np.any(below | above):
        raise ValueError("Every template age must fall inside age_bin_edges.")
    return ids


def smoothness_rows(ages: NDArray[np.float64]) -> NDArray[np.float64]:
    """First-difference rows, weighted so distant ages barely couple."""
    n = ages.size
    if n < 2:
        return np.zeros((0, n))
    order = np.argsort(ages, kind="mergesort")
    rows = []
    logt = np.log10(ages)
    for i in range(n - 1):
        j0 = int(order[i])
        j1 = int(order[i + 1])
        couple = float(np.exp(-abs(logt[j0] - logt[j1]) / _AGE_DEX_SCALE))
        if couple < _MIN_COUPLE:
            continue
        row = np.zeros(n)
        row[j0] = couple
        row[j1] = -couple
        rows.append(row)
    if not rows:
        return np.zeros((0, n))
    return np.vstack(rows)


@dataclass
class XRegularizer:
    mode: str
    ages: Optional[NDArray[np.float64]]
    strength: float = 1.0
    bin_edges: Optional[Tuple[float, ...]] = None

    @classmethod
    def from_config(
        cls,
        mode: Optional[str],
        template_ages: Optional[ArrayLike],
        n_comp: int,
        strength: float,
        bin_edges: Optional[Sequence[float]],
    ) -> Optional["XRegularizer"]:
        if mode is None:
            return None
        key = str(mode).lower()
        if key not in ("smooth_age", "age_bins"):
            raise ValueError("regularize_x must be None, 'smooth_age', or 'age_bins'.")
        if template_ages is None:
            raise ValueError("template_ages is required when regularize_x is set.")
        ages = _validate_ages(template_ages, n_comp)
        if key == "smooth_age" and strength <= 0:
            raise ValueError("regularize_strength must be positive for smooth_age.")
        edges = tuple(float(x) for x in bin_edges) if bin_edges is not None else None
        if key == "age_bins" and edges is not None:
            age_bin_ids(ages, edges)
        return cls(mode=key, ages=ages, strength=float(strength), bin_edges=edges)

    def restrict(self, keep: NDArray[np.bool_]) -> "XRegularizer":
        if self.ages is None:
            return self
        return XRegularizer(
            mode=self.mode,
            ages=self.ages[keep],
            strength=self.strength,
            bin_edges=self.bin_edges,
        )

    def solve(
        self,
        weighted_design: NDArray[np.float64],
        weighted_obs: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        n_comp = weighted_design.shape[1]
        if self.mode == "age_bins":
            return self._solve_bins(weighted_design, weighted_obs, n_comp)
        D = smoothness_rows(self.ages)
        if D.size == 0:
            x_hat, _ = nnls(weighted_design, weighted_obs)
            return x_hat
        used = weighted_obs != 0
        med_w = float(np.median(np.abs(weighted_obs[used]))) if np.any(used) else 1.0
        scale = self.strength * max(med_w, 1e-8)
        a_aug = np.vstack([weighted_design, scale * D])
        b_aug = np.concatenate([weighted_obs, np.zeros(D.shape[0])])
        x_hat, _ = nnls(a_aug, b_aug)
        return x_hat

    def _solve_bins(
        self,
        weighted_design: NDArray[np.float64],
        weighted_obs: NDArray[np.float64],
        n_comp: int,
    ) -> NDArray[np.float64]:
        ids = age_bin_ids(self.ages, self.bin_edges)
        bins = np.unique(ids)
        collapsed = np.empty((weighted_design.shape[0], bins.size))
        members = []
        for k, b in enumerate(bins):
            m = np.flatnonzero(ids == b)
            members.append(m)
            collapsed[:, k] = weighted_design[:, m].mean(axis=1)
        y, _ = nnls(collapsed, weighted_obs)
        x = np.zeros(n_comp)
        for yk, m in zip(y, members):
            x[m] = yk / m.size
        return x
