# Fiber Fundamental Mode Solver

This program computes the fundamental mode (LP01) of a step-index optical fiber by solving the Bessel functions of the fiber using a dichotomy method.  
It outputs the **effective index** (neff), the **mode field diameters (MFDs)**, the **effective area** (Aeff), the **V-number**, and optionally the **chromatic dispersion** and **group velocity**.

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
A_{\mathrm{eff}} \approx \pi \ w_0^2
$$

- **Rigorous field integral**

$$
A_{\mathrm{eff}} = \frac{\left( 2\pi \int_0^\infty |\psi(r)|^2 r \, dr \right)^2}
{2\pi \int_0^\infty |\psi(r)|^4 r \, dr}
$$

---

## Chromatic Dispersion and Group Velocity

The solver can also compute and plot the spectral dependence of:

- **Effective index** $n_\mathrm{eff}(\lambda)$  
- **Group velocity** $v_g(\lambda)$  
- **Chromatic dispersion** $D(\lambda)$ (in ps/nm/km)

The formulas used are:

- **Group index**:

$$
n_g(\lambda) = n_\mathrm{eff}(\lambda) - \lambda \, \frac{dn_\mathrm{eff}}{d\lambda}
$$

- **Group velocity**:

$$
v_g(\lambda) = \frac{c}{n_g(\lambda)}
$$

- **Chromatic dispersion**:

$$
D(\lambda) = -\frac{\lambda}{c} \, \frac{d^2 n_\mathrm{eff}}{d\lambda^2}
$$

where $c$ is the speed of light in vacuum.

⚠️ To reduce numerical noise, the first and last points (for group velocity) and the two first/last points (for dispersion) are discarded in the plots.

---

## V-Number

The solver computes the **V-number**:

$$
V = \frac{2\pi a}{\lambda} \, \mathrm{NA}
$$

where *a* is the fiber core radius and NA the numerical aperture.

---

## Advantages

Compared to simple approximations (e.g. the **Marcus approximation**, which can be inaccurate for practical fibers),  
this solver provides **higher precision** while remaining significantly **simpler to use** than a full vectorial electromagnetic mode solver.  
It also allows direct visualization of the **mode profile**, **MFDs**, and the **dispersion characteristics** of the fiber.

---

## References

- *Introduction to Fiber Optics*, Ajoy Ghatak and K. Thyagarajan  
- ITU-T G.650.1, *Definitions and test methods for mode field parameters*  
- G. P. Agrawal, *Nonlinear Fiber Optics* (for dispersion definitions)  


<img width="963" height="802" alt="image" src="https://github.com/user-attachments/assets/e26d0342-f6e2-4d78-92f9-d0c15b512c59" />
<img width="787" height="454" alt="image" src="https://github.com/user-attachments/assets/ef40b564-fdd0-434d-a501-a3e7c46f38c3" />

