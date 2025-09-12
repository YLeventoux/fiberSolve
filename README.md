# Fiber Fundamental Mode Solver

This program computes the fundamental mode (LP01) of a step-index optical fiber by solving the Bessel functions of the fiber using a dichotomy method.  
It outputs the **effective index** (neff), the **mode field diameters (MFDs)**, and the **effective area** (Aeff).

---

## Mode Field Diameter (MFD)

The mode field diameter characterizes the transverse extent of the fundamental mode and is often more relevant than the core diameter.  
Three definitions are implemented:

- **Gaussian approximation (1/e² criterion)**  
  Approximates the fundamental mode by a Gaussian and defines the MFD at the 1/e² intensity points.

- **Near-field rms definition**

$$
d_N = 2 \sqrt{ \frac{2 \int_0^\infty r^3 \, |\psi(r)|^2 \, dr}{\int_0^\infty r \, |\psi(r)|^2 \, dr} }
$$

where ψ(r) is the radial field profile.  
This is also referred to as the *4σ definition*.

- **Petermann II definition**

$$
w_P = \sqrt{ \frac{2 \int_0^\infty |\psi(r)|^2 r \, dr}{\int_0^\infty \left(\frac{d\psi}{dr}\right)^2 r \, dr} }
$$

$$
d_P = 2 \, w_P
$$

---

## Effective Area (Aeff)

Two approaches are available:

- **Gaussian approximation**

$$
A_{\mathrm{eff}} \approx \pi \, w_0^2
$$

- **Rigorous field integral**

$$
A_{\mathrm{eff}} = \frac{\left( 2\pi \int_0^\infty |\psi(r)|^2 r \, dr \right)^2}
{2\pi \int_0^\infty |\psi(r)|^4 r \, dr}
$$

---

## Advantages

Compared to simple approximations (e.g. the **Marcus approximation**, which can be inaccurate for practical fibers),  
this solver provides **higher precision** while remaining significantly **simpler to use** than a full vectorial electromagnetic mode solver.

---

## References

- *Introduction to Fiber Optics*, Ajoy Ghatak and K. Thyagarajan  
- ITU-T G.650.1, *Definitions and test methods for mode field parameters*




<img width="963" height="802" alt="image" src="https://github.com/user-attachments/assets/e26d0342-f6e2-4d78-92f9-d0c15b512c59" />
<img width="787" height="454" alt="image" src="https://github.com/user-attachments/assets/ef40b564-fdd0-434d-a501-a3e7c46f38c3" />

