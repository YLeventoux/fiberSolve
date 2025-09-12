# -*- coding: utf-8 -*-
"""
Fundamental mode solver of step-index fiber
with neff, MFDs, Aeff, V-number and dispersion calculation

Adapted from MATLAB code by Yann Leventoux (2025)
"""

import tkinter as tk
from tkinter import ttk, scrolledtext
import numpy as np
from scipy.special import j0, j1, k0 as besselk0, k1 as besselk1
from scipy.integrate import simps
import matplotlib.pyplot as plt

# -------- Sellmeier ----------
def silica_index(lambda_m, XGe=0):
    GeB1, GeB2, GeB3 = 0.80686642, 0.71815848, 0.85416831
    GeC1, GeC2, GeC3 = (0.068972606**2)*1e-12, (0.15396605**2)*1e-12, (11.841931**2)*1e-12
    SiB1, SiB2, SiB3 = 0.696166300, 0.407942600, 0.897479400
    SiC1, SiC2, SiC3 = 4.67914826e-3*1e-12, 1.35120631e-2*1e-12, 97.9340025*1e-12
    
    val = 1 + ((SiB1+XGe*(GeB1-SiB1))*lambda_m**2)/(lambda_m**2-(np.sqrt(SiC1)+XGe*(np.sqrt(GeC1)-np.sqrt(SiC1)))**2) \
            + ((SiB2+XGe*(GeB2-SiB2))*lambda_m**2)/(lambda_m**2-(np.sqrt(SiC2)+XGe*(np.sqrt(GeC2)-np.sqrt(SiC2)))**2) \
            + ((SiB3+XGe*(GeB3-SiB3))*lambda_m**2)/(lambda_m**2-(np.sqrt(SiC3)+XGe*(np.sqrt(GeC3)-np.sqrt(SiC3)))**2)
    return np.sqrt(val)

# -------- Equation de dispersion LP01 ----------
def eq_disp(n1,n2,rco,neff,lam):
    k0 = 2*np.pi/lam
    beta = k0*neff
    U = rco*np.sqrt((k0*n1)**2 - beta**2)
    W = rco*np.sqrt(beta**2 - (k0*n2)**2)
    lhs = j0(U)/(U*j1(U))
    rhs = besselk0(W)/(W*besselk1(W))   # corrigé
    return lhs-rhs

# -------- Dichotomie pour trouver neff ----------
def solve_neff(n1,n2,rco,lam,Nguess=100):
    neffVect = np.linspace(n2,n1,10000)
    Veff = 2*np.pi/lam*rco*np.sqrt(n1**2-neffVect**2)
    neffMin = neffVect[np.argmin(np.abs(Veff-2.405))]
    a,b = n1, neffMin
    for _ in range(Nguess):
        n1g = a+(b-a)/3
        n2g = b-(b-a)/3
        f1 = eq_disp(n1,n2,rco,n1g,lam)
        f2 = eq_disp(n1,n2,rco,n2g,lam)
        if f1*f2<0:
            b = n2g
        else:
            a = n1g
    return n1g

# -------- LP01 champ ----------
def besselLP01(r,n1,n2,rco,neff,lam):
    k0 = 2*np.pi/lam
    beta = k0*neff
    kT = np.sqrt(n1**2*k0**2-beta**2)
    gamma = np.sqrt(beta**2-n2**2*k0**2)
    coef = j0(kT*rco)/besselk0(gamma*rco)   # corrigé
    U = np.zeros_like(r)
    mask = np.abs(r)<rco
    U[mask] = j0(kT*r[mask])
    U[~mask] = coef*np.real(besselk0(gamma*np.abs(r[~mask])))  # corrigé
    return U

# -------- Mode solver ----------
def computeMode(a,rcl,n1,n2,lam):
    r = np.linspace(-rcl,rcl,3000)
    neff = solve_neff(n1,n2,a,lam)
    U = besselLP01(r,n1,n2,a,neff,lam)
    rp,Ep = r[r>=0], U[r>=0]
    normFactor = np.sqrt(simps(np.abs(Ep)**2*rp,rp)*2*np.pi)
    Ep = Ep/normFactor

    I = np.abs(Ep)**2/np.max(np.abs(Ep)**2)
    W0 = np.max(rp[I>0.1353])
    MFD_gauss = 2*W0

    dEp = np.gradient(Ep,rp)
    num = simps(rp*Ep**2,rp)
    den = simps(rp*dEp**2,rp)
    wp = np.sqrt(2*num/den)
    MFD_pet = 2*wp

    num = simps((rp**3)*np.abs(Ep)**2,rp)
    den = simps(rp*np.abs(Ep)**2,rp)
    wsig = 2*np.sqrt(num/den)
    MFD_4sigma = 2*wsig/np.sqrt(2)

    numA = (simps((np.abs(Ep)**2)*rp,rp))**2
    denA = simps((np.abs(Ep)**4)*rp,rp)
    Aeff = 2*np.pi*numA/denA

    return MFD_gauss,MFD_pet,MFD_4sigma,Aeff,neff,rp,Ep

# -------- Dispersion ----------
def calc_dispersion(paramType,paramValue,a):
    lamVec = np.linspace(0.55e-6,2.4e-6,100)
    neff_list = []
    for lam in lamVec:
        n2 = silica_index(lam,0)
        if paramType=="NA":
            n1 = np.sqrt(n2**2+paramValue**2)
        elif paramType=="delta n":
            n1 = n2+paramValue
        else:
            n1 = paramValue
        neff = solve_neff(n1,n2,a,lam)
        neff_list.append(neff)
    neffArr = np.array(neff_list)
    dneff = np.gradient(neffArr,lamVec)
    d2neff = np.gradient(dneff,lamVec)
    c=3e8
    ng = neffArr-lamVec*dneff
    vg = c/ng
    D = -(lamVec/c)*d2neff
    D_ps = D*1e6

    # Figures
    plt.figure(2); plt.clf()
    plt.plot(lamVec*1e6,neffArr,'-b')
    plt.xlabel("λ (µm)"); plt.ylabel("n_eff"); plt.title("n_eff vs λ"); plt.grid()

    plt.figure(3); plt.clf()
    plt.plot(lamVec[1:-1]*1e6,vg[1:-1]/1e8,'-g')  # enlève 1er et dernier point
    plt.xlabel("λ (µm)"); plt.ylabel("v_g (10^8 m/s)"); plt.title("Group velocity"); plt.grid()

    plt.figure(4); plt.clf()
    plt.plot(lamVec[2:-2]*1e6,D_ps[2:-2],'-r')  # enlève 2 premiers et 2 derniers points
    plt.xlabel("λ (µm)"); plt.ylabel("D (ps/nm/km)"); plt.title("Chromatic dispersion"); plt.grid()

    return lamVec,D_ps

# -------- GUI ----------
class FiberGUI:
    def __init__(self,root):
        self.root=root
        root.title("Fundamental mode calculator")

        tk.Label(root,text="Core diameter (µm):").grid(row=0,column=0,sticky="w")
        self.dCore=tk.DoubleVar(value=8.2); tk.Entry(root,textvariable=self.dCore).grid(row=0,column=1)

        tk.Label(root,text="Wavelength (µm):").grid(row=1,column=0,sticky="w")
        self.lam=tk.DoubleVar(value=1.55); tk.Entry(root,textvariable=self.lam).grid(row=1,column=1)

        tk.Label(root,text="Given parameter:").grid(row=2,column=0,sticky="w")
        self.typeMenu=ttk.Combobox(root,values=["NA","delta n","n1"]); self.typeMenu.current(0); self.typeMenu.grid(row=2,column=1)

        tk.Label(root,text="Value:").grid(row=3,column=0,sticky="w")
        self.val=tk.DoubleVar(value=0.14); tk.Entry(root,textvariable=self.val).grid(row=3,column=1)

        self.checkDisp=tk.BooleanVar(value=False)
        tk.Checkbutton(root,text="neff(λ), v_g, dispersion",variable=self.checkDisp).grid(row=4,column=0,columnspan=2)

        tk.Button(root,text="Calculate",command=self.calculate).grid(row=5,column=0,columnspan=2,pady=5)

        self.output=scrolledtext.ScrolledText(root,width=50,height=12)
        self.output.grid(row=6,column=0,columnspan=2)

    def calculate(self):
        a=self.dCore.get()/2*1e-6
        lam=self.lam.get()*1e-6
        n2=silica_index(lam,0)
        val=self.val.get()
        if self.typeMenu.get()=="NA": n1=np.sqrt(n2**2+val**2)
        elif self.typeMenu.get()=="delta n": n1=n2+val
        else: n1=val
        NA=np.sqrt(n1**2-n2**2)
        rcl=5*a
        MFDg,MFDp,MFD4,Aeff_rig,neff,rp,Ep=computeMode(a,rcl,n1,n2,lam)
        Aeff_apx=np.pi*(MFDg/2)**2
        V=2*np.pi*a*NA/lam

        self.output.delete("1.0",tk.END)
        self.output.insert(tk.END,
            f"n1={n1:.6f}, n2={n2:.6f}\n"
            f"neff={neff:.6f}\n"
            f"MFD (Gaussian 1/e²)={MFDg*1e6:.3f} µm\n"
            f"MFD (Petermann II) ={MFDp*1e6:.3f} µm\n"
            f"MFD (4σ rms)       ={MFD4*1e6:.3f} µm\n"
            f"Aeff (πw0²)        ={Aeff_apx*1e12:.3f} µm²\n"
            f"Aeff (rigorous)    ={Aeff_rig*1e12:.3f} µm²\n"
            f"V-number={V:.3f}\n")

        # ---- Plot mode ----
        plt.figure(1); plt.clf()
        plt.plot(rp*1e6,np.abs(Ep)**2/np.max(np.abs(Ep)**2),'b',label="Mode")
        plt.axvline(MFDg/2*1e6,color='r',ls='--',label="W0 Gaussian")
        plt.axvline(MFDp/2*1e6,color='g',ls='--',label="W0 Petermann")
        plt.axvline(MFD4/2*1e6,color='m',ls='--',label="W0 4σ")
        plt.xlabel("Radius (µm)"); plt.ylabel("Norm. intensity"); plt.title("LP01 mode"); plt.legend(); plt.grid()

        if self.checkDisp.get():
            lamVec,Dps=calc_dispersion(self.typeMenu.get(),val,a)
            # Interpolation au lambda choisi
            D_interp=np.interp(lam,lamVec,Dps)
            self.output.insert(tk.END,f"Dispersion at {lam*1e6:.2f} µm = {D_interp:.2f} ps/nm/km\n")

        plt.draw(); plt.pause(0.001)

# -------- Run ----------
if __name__=="__main__":
    root=tk.Tk(); app=FiberGUI(root); root.mainloop()
