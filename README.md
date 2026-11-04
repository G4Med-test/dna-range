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

- the **track length** of the primary;
- the **projected range** along the initial direction;
- the **penetration** (distance from the origin to the end point).

| Macro | Particle | Physics | Energies | Events per energy |
| --- | --- | --- | --- | --- |
| `e-_dna_opt2.mac` | e- | `dna_opt2` | 12 eV – 0.999 MeV | 10000 |
| `e-_dna_opt4.mac` | e- | `dna_opt4` | 12 eV – 9.999 MeV | 10000 |
| `e-_dna_opt6.mac` | e- | `dna_opt6` | 12 eV – 200 keV | 10000 |
| `e-_dna_opt8.mac` | e- | `dna_opt8` | 12 eV – 0.999 MeV | 10000 |
| `proton_dna_opt2.mac` | proton | `dna_opt2` | 100 eV – 300 MeV | 100 |
| `proton_dna_opt4.mac` | proton | `dna_opt4` | 100 eV – 300 MeV | 100 |
| `proton_dna_opt6.mac` | proton | `dna_opt6` | 100 eV – 300 MeV | 100 |
| `proton_dna_opt8.mac` | proton | `dna_opt8` | 100 eV – 300 MeV | 100 |
| `alpha_dna_opt2.mac` | alpha | `dna_opt2` | 100 eV – 400 MeV | 100 |
| `alpha_dna_opt4.mac` | alpha | `dna_opt4` | 100 eV – 400 MeV | 100 |
| `alpha_dna_opt6.mac` | alpha | `dna_opt6` | 100 eV – 400 MeV | 100 |
| `alpha_dna_opt8.mac` | alpha | `dna_opt8` | 100 eV – 400 MeV | 100 |

`validation/parser.py` publishes one plot per quantity and particle (x: energy in
MeV, y: length in nm, error: standard deviation / √N). The track length plot also
carries the **ICRU 90 CSDA range**.

## 📚 Reference: ICRU 90

`validation/icru90.py` returns the tabulated ICRU90 range values.

## Import notes

- `main.cc`, `include/`, `src/` are `range.cc`, `include/` and `src/` of the
  Geant4 v11.3.2 example (the version of the default base image), unchanged.
- `macro/range.in` and `macro/histo.in` are the example's macros, kept for
  reference; the validation macros are the `*_dna_opt*.mac` files listed in
  `validation/config.json`. `docs/geant4-example-README.txt` is the example's
  README. `plot.C` is replaced by `validation/parser.py`.
- `macro/unit.mac` runs one 1 keV electron on every push.

For the CI, the validation flow and running everything locally, see the
[template README](https://github.com/G4Med-test/template#readme).

Geant4-DNA references for this example: Med. Phys. 51 (2024) 5873–5889;
Med. Phys. 45 (2018) e722–e739; Phys. Med. 31 (2015) 861–874; Med. Phys. 37
(2010) 4692–4708; Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178.
