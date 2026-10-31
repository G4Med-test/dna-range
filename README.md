# dna-range

**Range of electrons, protons and alpha particles in liquid water** with the
Geant4-DNA track-structure physics, compared with the CSDA ranges derived from
ICRU Report 90. Built on the Geant4 example
[`extended/medical/dna/range`](https://github.com/Geant4/geant4/tree/v11.3.2/examples/extended/medical/dna/range).

---

## 🔬 What this test does

Mono-energetic particles start at the centre of a 1 m radius sphere of
`G4_WATER`, in a random direction. Secondary electrons are killed at their first
step, so only the primary is tracked, down to the 10 eV tracking cut. For each
energy the application writes to `range.txt` the mean and standard deviation over
the events of:

- the **track length** of the primary (for protons and alphas, summed over the
  Geant4-DNA charge states: H, He+, He);
- the **projected range** along the initial direction;
- the **penetration** (distance from the origin to the end point).

| Macro | Particle | Physics | Energies | Events per energy |
| --- | --- | --- | --- | --- |
| `e-_dna_opt2.mac` | e- | `dna_opt2` | 12 eV – 1 MeV | 1000 |
| `proton_dna_opt2.mac` | proton | `dna_opt2` | 1 keV – 10 MeV | 1000 (200 at 5, 10 MeV) |
| `alpha_dna_opt2.mac` | alpha | `dna_opt2` | 1 keV – 10 MeV | 1000 (200 at 5, 10 MeV) |

`validation/parser.py` publishes one plot per quantity and particle (x: energy in
MeV, y: length in nm, error: standard deviation / √N). The track length plot also
carries the **ICRU 90 CSDA range**, from 10 keV (e-), 100 keV (protons) and
500 keV (alphas) up, where it does not depend on the low-energy assumptions
below. Projected range and penetration have no reference and are compared between
Geant4 versions.

## 📚 Reference: ICRU 90

`validation/icru90.py` (identical in `dna-spower`) computes the CSDA range in
liquid water by integrating the reciprocal of the ICRU 90 stopping power from
1 keV, plus an estimate of the range at 1 keV:

- **protons, alphas**: electronic stopping power of `G4_WATER` from
  `G4ICRU90StoppingData` (the ICRU 90 table embedded in Geant4,
  `source/materials/src/G4ICRU90StoppingData.cc`), log-log interpolated. Nuclear
  stopping is not included, hence the higher lower limits; below 1 keV the
  stopping power is taken proportional to the velocity.
- **electrons**: collision stopping power from the Bethe formula (ICRU 37
  eq. 2.1) with the ICRU 90 mean excitation energy I = 78 eV and density-effect
  correction (zero below 0.5 MeV), as ESTAR computes it (with I = 75 eV it
  reproduces ESTAR to 4 digits). Geant4-DNA has no bremsstrahlung, so the
  collision, not the total, stopping power is the right reference. Below 1 keV
  the stopping power is taken as constant.

The track length is the path length the CSDA range approximates; it includes no
detour factor.

## Import notes

- `main.cc`, `include/`, `src/` are `range.cc`, `include/` and `src/` of the
  Geant4 v11.3.2 example (the version of the default base image), unchanged.
- `macro/range.in` and `macro/histo.in` are the example's macros, kept for
  reference; the validation macros are the `*_dna_opt2.mac` files listed in
  `validation/config.json`. `docs/geant4-example-README.txt` is the example's
  README. `plot.C` is replaced by `validation/parser.py`.
- `macro/unit.mac` runs one 1 keV electron on every push.

For the CI, the validation flow and running everything locally, see the
[template README](https://github.com/G4Med-test/template#readme).

Geant4-DNA references for this example: Med. Phys. 51 (2024) 5873–5889;
Med. Phys. 45 (2018) e722–e739; Phys. Med. 31 (2015) 861–874; Med. Phys. 37
(2010) 4692–4708; Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178.
