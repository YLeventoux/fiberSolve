# -*- coding: utf-8 -*-
"""
solving of fundamental mode of a given fiber (by solving bessel function)

Created on Wed Sep 10 16:03:56 2025

@author: Leventoux Yann
"""

import tkinter as tk
from tkinter import ttk, scrolledtext
import numpy as np
from scipy.special import j0, k0 as besselk0
from scipy.integrate import simps
import matplotlib.pyplot as plt

# -------- Sellmeier for silica / Ge-doped silica ----------
def silica_index(lambda_m, XGe=0):
    GeB1, GeB2, GeB3 = 0.80686642, 0.71815848, 0.85416831
    GeC1, GeC2, GeC3 = (0.068972606**2)*1e-12, (0.15396605**2)*1e-12, (11.841931**2)*1e-12
    SiB1, SiB2, SiB3 = 0.696166300, 0.407942600, 0.897479400
    SiC1, SiC2, SiC3 = 4.67914826e-3*1e-12, 1.35120631e-2*1e-12, 97.9340025*1e-12
    
    val = 1 + ((SiB1+XGe*(GeB1-SiB1))*lambda_m**2)/(lambda_m**2-(np.sqrt(SiC1)+XGe*(np.sqrt(GeC1)-np.sqrt(SiC1)))**2) \
            + ((SiB2+XGe*(GeB2-SiB2))*lambda_m**2)/(lambda_m**2-(np.sqrt(SiC2)+XGe*(np.sqrt(GeC2)-np.sqrt(SiC2)))**2) \
            + ((SiB3+XGe*(GeB3-SiB3))*lambda_m**2)/(lambda_m**2-(np.sqrt(SiC3)+XGe*(np.sqrt(GeC3)-np.sqrt(SiC3)))**2)
    return np.sqrt(val)

# -------- LP01 field -----------
def besselLP01(r, n1, n2, rco, neff, wl):
    k0 = 2*np.pi/wl
    beta = 2*np.pi*neff/wl
    kT = np.sqrt(n1**2*k0**2 - beta**2)
    gamma = np.sqrt(beta**2 - n2**2*k0**2)
    coef = j0(kT*rco)/besselk0(gamma*rco)
    U = np.zeros_like(r)
    mask_core = np.abs(r) < rco
    U[mask_core] = j0(kT*r[mask_core])
    U[~mask_core] = coef*np.real(besselk0(gamma*np.abs(r[~mask_core])))
    return U

# -------- Mode solver -----------
def computeMode(a, rcl, n1, n2, lam):
    r = np.linspace(-rcl, rcl, 3000)
    neffVect = np.linspace(n2, n1, 10000)
    Veff = 2*np.pi/lam*a*np.sqrt(n1**2 - neffVect**2)
    neffMin = neffVect[np.argmin(np.abs(Veff-2.405))]
    
    a1, b1 = n1, neffMin
    for _ in range(100):
        nGuess1 = a1+(b1-a1)/3
        nGuess2 = b1-(b1-a1)/3
        U1 = besselLP01(r,n1,n2,a,nGuess1,lam)
        U2 = besselLP01(r,n1,n2,a,nGuess2,lam)
        U1p = np.gradient(np.gradient(U1, r), r)
        U2p = np.gradient(np.gradient(U2, r), r)
        if np.min(U1p) > np.min(U2p):
            b1 = nGuess2
        else:
            a1 = nGuess1
    neff = nGuess1

    U = besselLP01(r,n1,n2,a,neff,lam)
    rp = r[r>=0]
    Ep = U[r>=0]

    normFactor = np.sqrt(simps(np.abs(Ep)**2 * rp, rp) * 2*np.pi)
    Ep = Ep/normFactor

    I = np.abs(Ep)**2
    I = I/np.max(I)
    temp = rp[I>0.1353]
    W0 = np.max(temp)
    MFD_gauss = 2*W0

    num = simps((rp**2)*np.abs(Ep)**2, rp)
    num = num**2
    den = simps((rp**3)*np.abs(Ep)**4, rp)
    wp = np.sqrt(2*num/den)
    MFD_pet = 2*wp

    num2 = simps(np.abs(Ep)**2*rp, rp)
    num2 = num2**2
    den2 = simps(np.abs(Ep)**4*rp, rp)
    Aeff = 2*np.pi*num2/den2

    return MFD_gauss, MFD_pet, Aeff, neff, rp, Ep

# -------- GUI ----------
class FiberGUI:
    def __init__(self, root):
        self.root = root
        root.title("Fundamental mode calculator")

        tk.Label(root,text="Core diameter (µm):").grid(row=0,column=0,sticky="w")
        self.dCore = tk.DoubleVar(value=8.2)
        tk.Entry(root,textvariable=self.dCore).grid(row=0,column=1)

        tk.Label(root,text="Wavelength (µm):").grid(row=1,column=0,sticky="w")
        self.lambdaVar = tk.DoubleVar(value=1.55)
        tk.Entry(root,textvariable=self.lambdaVar).grid(row=1,column=1)

        tk.Label(root,text="Given parameter:").grid(row=2,column=0,sticky="w")
        self.typeMenu = ttk.Combobox(root,values=["NA","delta n","n1"])
        self.typeMenu.current(0)
        self.typeMenu.grid(row=2,column=1)

        tk.Label(root,text="Value:").grid(row=3,column=0,sticky="w")
        self.valueVar = tk.DoubleVar(value=0.14)
        tk.Entry(root,textvariable=self.valueVar).grid(row=3,column=1)

        tk.Button(root,text="Calculate",command=self.calculate).grid(row=4,column=0,columnspan=2,pady=5)

        self.output = scrolledtext.ScrolledText(root,width=50,height=10)
        self.output.grid(row=5,column=0,columnspan=2)

    def calculate(self):
        a = self.dCore.get()/2*1e-6
        lam = self.lambdaVar.get()*1e-6
        n2 = silica_index(lam,0)
        val = self.valueVar.get()

        if self.typeMenu.get()=="NA":
            NA = val
            n1 = np.sqrt(n2**2+NA**2)
        elif self.typeMenu.get()=="delta n":
            n1 = n2+val
        elif self.typeMenu.get()=="n1":
            n1 = val
        else:
            n1 = n2+0.001

        NA = np.sqrt(n1**2-n2**2)
        rcl = 10*a

        MFD_gauss, MFD_pet, Aeff_rig, neff, rp, Ep = computeMode(a,rcl,n1,n2,lam)
        Aeff_approx = np.pi*(MFD_gauss/2)**2
        V = 2*np.pi*a*NA/lam

        self.output.delete("1.0",tk.END)
        self.output.insert(tk.END,
            f"n1 = {n1:.6f}, n2 = {n2:.6f}\n"
            f"neff = {neff:.6f}\n"
            f"MFD (Gaussian 1/e²) = {MFD_gauss*1e6:.3f} µm\n"
            f"MFD (Petermann II)  = {MFD_pet*1e6:.3f} µm\n"
            f"Aeff (pi*w0²)       = {Aeff_approx*1e12:.3f} µm²\n"
            f"Aeff (field int.)   = {Aeff_rig*1e12:.3f} µm²\n"
            f"V-number = {V:.3f}\n")

        # ---- Separate plot (with refresh) ----
        plt.figure(1, figsize=(6,4))
        plt.clf()
        plt.plot(rp*1e6, np.abs(Ep)**2/np.max(np.abs(Ep)**2),
                 'b', label="Mode intensity")
        plt.axvline(MFD_gauss/2*1e6, color='r', linestyle='--', label="W0 Gaussian")
        plt.axvline(MFD_pet/2*1e6, color='g', linestyle='--', label="W0 Petermann")
        plt.xlabel("Radius (µm)")
        plt.ylabel("Normalized intensity")
        plt.title("LP01 mode profile and MFDs")
        plt.legend()
        plt.grid()

        plt.draw()
        plt.pause(0.001)   # <--- indispensable pour rafraîchir

# -------- Run ----------
if __name__=="__main__":
    root = tk.Tk()
    app = FiberGUI(root)
    root.mainloop()


