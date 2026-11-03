"""Finite circular-source / square-detector quadrature and input validation."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import least_squares, brentq
from scipy.special import xlogy
from scipy.stats import chi2
from numpy.polynomial.legendre import leggauss

RS = .004875
SIDE = .010
DELTA0, DELTA_SD = .018, .0015

def efficiency(r, radial=24, angular=96, derivative=False):
    r = np.atleast_1d(r).astype(float)
    if np.any(r <= 0):
        raise ValueError('Source–detector separation must be positive')
    z, w = leggauss(radial)
    rho = RS * np.sqrt((z + 1) / 2)
    phi = 2*np.pi*(np.arange(angular)+.5)/angular
    xs = rho[:, None]*np.cos(phi)[None, :]
    ys = rho[:, None]*np.sin(phi)[None, :]
    x1, x2 = -SIDE/2-xs, SIDE/2-xs
    y1, y2 = -SIDE/2-ys, SIDE/2-ys
    rr = r[:, None, None]
    def corner(a, b):
        s = np.sqrt(rr*rr+a*a+b*b)
        if derivative:
            return -a*b*(2*rr*rr+a*a+b*b)/(s*(rr*rr+a*a)*(rr*rr+b*b))
        return np.arctan2(a*b, rr*s)
    omega = corner(x2,y2)-corner(x2,y1)-corner(x1,y2)+corner(x1,y1)
    return np.sum(omega * (w[None,:,None]/2), axis=(1,2))/(angular*4*np.pi)

def verify_geometry():
    r = np.array([.0265, .05, .118, .818])
    e = efficiency(r)
    fine = efficiency(r, 48, 192)
    h = 1e-7
    deriv = efficiency(r, derivative=True)
    fd = (efficiency(r+h)-efficiency(r-h))/(2*h)
    # Independent quadrature over the detector plane as well as source disc.
    z, w = leggauss(24); phi = 2*np.pi*(np.arange(96)+.5)/96
    rho = RS*np.sqrt((z+1)/2)
    sx = (rho[:,None]*np.cos(phi)).ravel()
    sy = (rho[:,None]*np.sin(phi)).ravel()
    ws = np.repeat(w/2/96, 96)
    q, v = leggauss(16)
    dx, dy = np.meshgrid(q*SIDE/2, q*SIDE/2, indexing='ij')
    wd = np.outer(v,v).ravel()*(SIDE/2)**2
    dsq = (dx.ravel()[None,:]-sx[:,None])**2+(dy.ravel()[None,:]-sy[:,None])**2
    direct = np.array([np.sum(ws[:,None]*wd[None,:]*rr/(rr*rr+dsq)**1.5)/(4*np.pi) for rr in r])
    checks = {'quadrature_max_relative_difference': float(np.max(abs(e/fine-1))),
              'detector_integral_max_relative_difference': float(np.max(abs(e/direct-1))),
              'derivative_max_relative_difference': float(np.max(abs(deriv/fd-1)))}
    if max(checks.values()) > 2e-6: raise AssertionError(checks)
    return checks

def read_data(path):
    df = pd.read_csv(path)
    for col in ['d','d_unc','n','t']:
        if col not in df: raise ValueError('Missing required column '+col)
    background = df['d'].astype(str).str.strip().str.lower().eq('background')
    omitted = df.loc[background].to_dict('records')
    data = df.loc[~background].copy()
    for c in ['d','d_unc','n','t']:
        data[c] = pd.to_numeric(data[c], errors='raise')
    if not np.isfinite(data[['d','d_unc','n','t']].to_numpy()).all(): raise ValueError('Nonfinite input')
    if ((data.d <= 0) | (data.d > 2) | (data.d_unc <= 0) | (data.d_unc >= data.d) | (data.n < 0) | (data.t <= 0)).any():
        raise ValueError('Implausible SI values or uncertainties; inspect input units')
    if not np.allclose(data.n, np.round(data.n)): raise ValueError('n must be integer event counts')
    return data.reset_index(drop=True), omitted
