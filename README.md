This program computes the fundamental mode (LP01) of a step-index optical fiber by solving the Bessel functions of the fiber using a dichotomy method.
It outputs the effective index (neff), the mode field diameter (waist), and the effective area (Aeff).

The mode field diameter can be evaluated in two ways:
• Gaussian approximation (1/e² criterion)
• Petermann II definition

The effective area (Aeff) can also be obtained with two approaches:
• Gaussian approximation (π·w₀²)
• Rigorous field integral formula

Compared to the Marcus approximation, which is often inaccurate for practical fibers, this method provides higher precision while remaining much simpler to use than a full electromagnetic mode solver.
<img width="963" height="802" alt="image" src="https://github.com/user-attachments/assets/e26d0342-f6e2-4d78-92f9-d0c15b512c59" />
