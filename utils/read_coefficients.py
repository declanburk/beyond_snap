import matplotlib.pyplot as plt
from plotting_tools import *

n = 10 # Order of polynomials
m = 0.01 # Evaluation discretisation

## Minimum snap plotting

P, T, B, k, J = read_spline_from_csv('min_snap_coefficients.csv')

print("Evaluating polynomials.")

Tt, Z = eval_normalised_poly_coeffs(P, T, k, n, m)

print("Plotting.")

font_size = 15
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': font_size})

plt.figure(1)
plt.plot(Tt, Z, c='k', ls='-', lw=2, label="Random Fixed-Time Trajectory")
plt.scatter(T, B, c="r")
plt.xlabel("Time (s)", fontsize=font_size)
plt.ylabel("Position (m)", fontsize=font_size)
plt.title("Snap %.2E" % J, fontsize=2*font_size)
plt.legend()

plt.show()