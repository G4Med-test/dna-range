"""ICRU 90 tabulated CSDA ranges in liquid water.

All tabulated values are embedded below; no external data files are needed.
Energy: MeV. Range: g/cm² (mass thickness).
Between tabulated energies, use log-log interpolation; no extrapolation.
"""

from bisect import bisect_left
from math import exp, log

# (kinetic energy in MeV, range in g/cm²), taken verbatim from supplied files.
RANGE_DATA = {
    'e-': (
        (0.001, 4.235E-6),
        (0.0015, 9.118E-6),
        (0.002, 1.524E-5),
        (0.003, 3.098E-5),
        (0.004, 5.110E-5),
        (0.005, 7.536E-5),
        (0.006, 1.036E-4),
        (0.008, 1.714E-4),
        (0.01, 2.537E-4),
        (0.015, 5.189E-4),
        (0.02, 8.632E-4),
        (0.03, 1.769E-3),
        (0.04, 2.939E-3),
        (0.05, 4.348E-3),
        (0.06, 5.978E-3),
        (0.08, 9.833E-3),
        (0.1, 1.439E-2),
        (0.15, 2.833E-2),
        (0.2, 4.512E-2),
        (0.3, 8.464E-2),
        (0.4, 1.294E-1),
        (0.5, 1.774E-1),
        (0.6, 2.275E-1),
        (0.8, 3.315E-1),
        (1, 4.384E-1),
        (1.5, 7.096E-1),
        (2, 9.811E-1),
        (3, 1.517E+00),
        (4, 2.041E+00),
        (5, 2.554E+00),
        (6, 3.057E+00),
        (8, 4.035E+00),
        (10, 4.980E+00),
        (15, 7.226E+00),
        (20, 9.327E+00),
        (30, 1.318E+01),
        (40, 1.667E+01),
        (50, 1.985E+01),
        (60, 2.277E+01),
        (80, 2.801E+01),
        (100, 3.259E+01),
        (150, 4.204E+01),
        (200, 4.955E+01),
        (300, 6.110E+01),
        (400, 6.986E+01),
        (500, 7.692E+01),
        (600, 8.284E+01),
        (800, 9.239E+01),
        (1000, 9.994E+01),
    ),
    'proton': (
        (0.0010, 6.319E-06),
        (0.0015, 8.969E-06),
        (0.0020, 1.137E-05),
        (0.0030, 1.560E-05),
        (0.0040, 1.930E-05),
        (0.0050, 2.262E-05),
        (0.0060, 2.567E-05),
        (0.0080, 3.113E-05),
        (0.0100, 3.599E-05),
        (0.0150, 4.657E-05),
        (0.0200, 5.578E-05),
        (0.0300, 7.187E-05),
        (0.0400, 8.613E-05),
        (0.0500, 9.935E-05),
        (0.0600, 1.120E-04),
        (0.0800, 1.364E-04),
        (0.1000, 1.607E-04),
        (0.1500, 2.249E-04),
        (0.2000, 2.967E-04),
        (0.3000, 4.645E-04),
        (0.4000, 6.640E-04),
        (0.5000, 8.945E-04),
        (0.6000, 1.155E-03),
        (0.8000, 1.765E-03),
        (1.0000, 2.487E-03),
        (1.5000, 4.753E-03),
        (2.0000, 7.639E-03),
        (3.0000, 1.514E-02),
        (4.0000, 2.482E-02),
        (5.0000, 3.656E-02),
        (6.0000, 5.028E-02),
        (8.0000, 8.346E-02),
        (10.0000, 1.240E-01),
        (15.0000, 2.558E-01),
        (20.0000, 4.289E-01),
        (30.0000, 8.910E-01),
        (40.0000, 1.498E+00),
        (50.0000, 2.240E+00),
        (60.0000, 3.111E+00),
        (80.0000, 5.212E+00),
        (100.0000, 7.759E+00),
        (150.0000, 1.586E+01),
        (200.0000, 2.609E+01),
        (300.0000, 5.170E+01),
        (400.0000, 8.263E+01),
        (500.0000, 1.176E+02),
        (600.0000, 1.556E+02),
        (800.0000, 2.382E+02),
        (1000.0000, 3.268E+02),
        (1500.0000, 5.625E+02),
        (2000.0000, 8.079E+02),
        (3000.0000, 1.307E+03),
        (4000.0000, 1.806E+03),
        (5000.0000, 2.302E+03),
        (6000.0000, 2.792E+03),
        (8000.0000, 3.758E+03),
        (10000.0000, 4.707E+03),
    ),
    'alpha': (
        (0.0010, 3.273E-06),
        (0.0015, 4.789E-06),
        (0.0020, 6.294E-06),
        (0.0030, 9.247E-06),
        (0.0040, 1.211E-05),
        (0.0050, 1.489E-05),
        (0.0060, 1.757E-05),
        (0.0080, 2.267E-05),
        (0.0100, 2.746E-05),
        (0.0150, 3.831E-05),
        (0.0200, 4.793E-05),
        (0.0300, 6.462E-05),
        (0.0400, 7.900E-05),
        (0.0500, 9.179E-05),
        (0.0600, 1.034E-04),
        (0.0800, 1.241E-04),
        (0.1000, 1.425E-04),
        (0.1500, 1.817E-04),
        (0.2000, 2.151E-04),
        (0.3000, 2.725E-04),
        (0.4000, 3.230E-04),
        (0.5000, 3.699E-04),
        (0.6000, 4.150E-04),
        (0.8000, 5.034E-04),
        (1.0000, 5.931E-04),
        (1.5000, 8.385E-04),
        (2.0000, 1.128E-03),
        (3.0000, 1.844E-03),
        (4.0000, 2.737E-03),
        (5.0000, 3.798E-03),
        (6.0000, 5.020E-03),
        (8.0000, 7.924E-03),
        (10.0000, 1.142E-02),
        (15.0000, 2.258E-02),
        (20.0000, 3.702E-02),
        (30.0000, 7.514E-02),
        (40.0000, 1.250E-01),
        (50.0000, 1.859E-01),
        (60.0000, 2.575E-01),
        (80.0000, 4.317E-01),
        (100.0000, 6.452E-01),
        (150.0000, 1.341E+00),
        (200.0000, 2.253E+00),
        (300.0000, 4.670E+00),
        (400.0000, 7.802E+00),
        (500.0000, 1.158E+01),
        (600.0000, 1.594E+01),
        (800.0000, 2.622E+01),
        (1000.0000, 3.830E+01),
    ),
}

PARTICLES = tuple(RANGE_DATA)

# Minimum kinetic energies (MeV) used by parser.py for reference curves.
# Preserve the original icru90.py validation thresholds.
RANGE_MIN_MEV = {"e-": 0.001, "proton": 0.001, "alpha": 0.001}


def range_table(particle):
    """Return (energies_MeV, ranges_g_cm2) as tuples."""
    try:
        rows = RANGE_DATA[particle]
    except KeyError:
        raise ValueError(f"Unknown particle {particle!r}; choose from {PARTICLES}") from None
    return tuple(e for e, _ in rows), tuple(r for _, r in rows)


def energy_limits(particle):
    """Return minimum and maximum tabulated energy (MeV)."""
    energies, _ = range_table(particle)
    return energies[0], energies[-1]


def covers(particle, energy):
    """Return True when energy is within tabulated limits."""
    low, high = energy_limits(particle)
    return low <= energy <= high


def csda_range(particle, energy):
    """Return CSDA mass range (g/cm²) for kinetic energy (MeV)."""
    energies, ranges = range_table(particle)
    if not energies[0] <= energy <= energies[-1]:
        raise ValueError(
            f"{energy:g} MeV is outside the table for {particle} "
            f"[{energies[0]:g}, {energies[-1]:g}] MeV"
        )

    i = bisect_left(energies, energy)
    if energies[i] == energy:
        return ranges[i]

    e0, e1 = energies[i - 1], energies[i]
    r0, r1 = ranges[i - 1], ranges[i]
    fraction = log(energy / e0) / log(e1 / e0)
    return exp(log(r0) + fraction * log(r1 / r0))


get_range = csda_range


if __name__ == "__main__":
    for particle in PARTICLES:
        print(f"{particle:7s} 1 MeV: {csda_range(particle, 1.0):.6g} g/cm²")
