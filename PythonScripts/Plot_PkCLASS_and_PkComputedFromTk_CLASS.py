#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

# ---------------- Top-hat window in Fourier space ----------------
def W_tophat(kR):
    kR = np.where(kR == 0, 1e-10, kR)  # avoid divide by zero
    return 3 * (np.sin(kR) - kR * np.cos(kR)) / kR**3

def sigma_r(k, Pk, R):
    """Compute sigma_R from 1D P(k)"""
    W = W_tophat(k * R)
    integrand = k**2 * Pk * W**2
    return np.sqrt(integrate.simpson(integrand, k) / (2*np.pi**2))

def compute_pk_from_Tk(tk_file, n_s, sigma8_target, h):
    """Compute P(k) from CLASS Tk file and normalize to sigma8"""
    data = np.loadtxt(tk_file, skiprows=8)
    k_hmpc = data[:, 0]           # k in h/Mpc
    T_tot = np.abs(data[:, 6])    # total matter transfer function

    # Unnormalized P(k)
    Pk_unnorm = T_tot**2 * k_hmpc**n_s

    # Normalize to sigma8
    s8_current = sigma_r(k_hmpc, Pk_unnorm, R=8.0)
    norm = sigma8_target**2 / s8_current**2
    Pk = norm * Pk_unnorm

    s8_final = sigma_r(k_hmpc, Pk, R=8.0)
    print(f"sigma_8 after normalization: {s8_final:.5f}, target: {sigma8_target}")
    print(f"Normalization factor applied: {norm:.5e}")

    return k_hmpc, Pk

def read_class_pk(class_file):
    """Read CLASS P(k) file"""
    data = np.loadtxt(class_file)
    k = data[:,0]
    Pk = data[:,1]
    return k, Pk

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file-class-Pk", type=str, required=True, help="CLASS P(k) file")
    parser.add_argument("--file-class-Tk", type=str, required=True, help="CLASS Tk.dat file")
    parser.add_argument("--n_s", type=float, required=True, help="Spectral index n_s")
    parser.add_argument("--sigma8", type=float, required=True, help="Target sigma8")
    parser.add_argument("--h", type=float, required=True, help="Hubble parameter h")
    parser.add_argument("--output", type=str, default="Pk_combined.png", help="Output plot filename")
    args = parser.parse_args()

    # Read CLASS P(k)
    k_class, Pk_class = read_class_pk(args.file_class_Pk)

    # Compute P(k) from Tk and normalize
    k_tk, Pk_tk = compute_pk_from_Tk(args.file_class_Tk, args.n_s, args.sigma8, args.h)

    # -------- Plot both spectra on the same figure --------
    plt.figure(figsize=(6,4))
    plt.loglog(k_class, Pk_class, 'ok', label="CLASS P(k)", markersize=3)
    plt.loglog(k_tk, Pk_tk, 'b-', label="Tk -> P(k), sigma8-normalized")
    plt.xlabel("$k$ [h/Mpc]")
    plt.ylabel("$P(k)$ [(Mpc/h)$^3$]")
    plt.legend()
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Saved combined plot to {args.output}")

if __name__ == "__main__":
    main()
