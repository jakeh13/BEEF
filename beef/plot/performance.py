# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Georgia Tech Research Corporation
"""Performance plots: insertion loss and splitting ratio."""
from __future__ import annotations

import os
from typing import List, Optional

import numpy as np
from matplotlib import pyplot as plt

from ..results.smatrix import SMatrix


def plot_performance(
    smat: SMatrix,
    output_dir: str,
    title: str = "",
    ports: Optional[List[str]] = None,
    filename: str = "performance.png",
) -> None:
    """Plot transmission vs wavelength for all ports.

    Args:
        smat: SMatrix result from any solver.
        output_dir: directory to save the PNG.
        title: plot title (e.g. device dimensions).
        ports: list of port names to plot; defaults to all ports in smat.
        filename: output filename.
    """
    if ports is None:
        ports = smat.ports()

    wls = smat.wavelengths
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    ax_lin, ax_db = axes

    for port in ports:
        T = smat.transmission(port=port)
        IL = smat.insertion_loss_db(port=port)
        ax_lin.plot(wls, T, label=port.capitalize())
        ax_db.plot(wls, IL, label=port.capitalize())

    ax_lin.set_xlabel("Wavelength (µm)")
    ax_lin.set_ylabel("Transmission")
    ax_lin.set_title(title or "Transmission")
    ax_lin.legend()
    ax_lin.grid(True)

    ax_db.set_xlabel("Wavelength (µm)")
    ax_db.set_ylabel("Insertion Loss (dB)")
    ax_db.set_title(title or "Insertion Loss")
    ax_db.legend()
    ax_db.grid(True)

    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, filename), dpi=300)
    plt.close("all")


def plot_splitting_ratio(
    smat: SMatrix,
    port_a: str,
    port_b: str,
    output_dir: str,
    title: str = "",
    filename: str = "splitting_ratio.png",
) -> None:
    """Plot power splitting ratio T_a / (T_a + T_b) vs wavelength.

    Args:
        smat: SMatrix result.
        port_a: numerator port label.
        port_b: denominator port label.
        output_dir: directory to save the PNG.
    """
    wls = smat.wavelengths
    ratio = smat.splitting_ratio(port_a, port_b)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(wls, ratio)
    ax.axhline(0.5, color="k", linestyle="--", label="Target 50/50")
    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel(f"Splitting Ratio ({port_a} / ({port_a}+{port_b}))")
    ax.set_title(title or "Splitting Ratio")
    ax.legend()
    ax.grid(True)
    ax.set_ylim(0, 1)
    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, filename), dpi=300)
    plt.close("all")
