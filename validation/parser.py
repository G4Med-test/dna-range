#!/usr/bin/env python

# dna-range portal parser: metadata() reads the macro, parse() the range.txt
# table that Run::EndOfRun appends after every /run/beamOn (one macro = one
# particle, one energy per run): energy (eV), then mean and standard deviation
# over the events (nm) of the track length, the projected range and the
# penetration of the primary particle.
#
# The track length of the primary (summed over its Geant4-DNA charge states for
# protons and alphas, see TrackingAction) is compared with the ICRU 90 CSDA range
# of liquid water (validation/icru90.py), from the energy where that range is
# reliable (icru90.RANGE_MIN_MEV) up. Projected range and penetration have no
# reference and are published for comparison between Geant4 versions.
import math
import os

import icru90
from geantval import energy_mev, getJSON, one_command

TEST = "dna-range"
PHYSLISTS = {"dna", "dna_opt1", "dna_opt2", "dna_opt3", "dna_opt4", "dna_opt5", "dna_opt6",
             "dna_opt7", "std_opt4"}
DENSITY = 1.0  # g/cm3, G4_WATER

# range.txt columns after the energy: (title, mean column, standard deviation column).
QUANTITIES = (("track length", 1, 2),
              ("projected range", 3, 4),
              ("penetration", 5, 6))


def energy_scan(commands):
    """Energy (MeV) and number of events of every /run/beamOn, in macro order."""
    energies, events = [], []
    energy = None
    for command, value in commands:
        if command == "/gun/energy":
            energy = energy_mev(value)
        elif command == "/run/beamOn":
            n = int(value)
            if energy is None or n < 2:
                raise ValueError("Each /run/beamOn needs an explicit /gun/energy before it "
                                 "and at least 2 events")
            energies.append(energy)
            events.append(n)
    if not energies:
        raise ValueError("No /run/beamOn in the macro")
    return energies, events


def metadata(commands):
    physlist = one_command(commands, "/range/phys/addPhysics")
    if physlist not in PHYSLISTS:
        # The application silently keeps its default for an unknown name.
        raise ValueError("Unknown physics constructor %r (expected one of %s)"
                         % (physlist, ", ".join(sorted(PHYSLISTS))))
    if one_command(commands, "/range/det/setMat") != "G4_WATER":
        raise ValueError("The ICRU 90 comparison needs /range/det/setMat G4_WATER")
    particle = one_command(commands, "/gun/particle")
    if particle not in icru90.PARTICLES:
        raise ValueError("Unsupported particle %r (expected one of %s)"
                         % (particle, ", ".join(icru90.PARTICLES)))
    energies, events = energy_scan(commands)
    return {"TEST": TEST, "PHYSLIST": physlist, "PARTICLE": particle,
            "TRACKING_CUT": one_command(commands, "/range/det/setTrackingCut"),
            "ENERGIES": energies, "NEVENTS": events}


def read_table(path, expected_energies):
    """Rows of range.txt, checked against the energies of the macro."""
    rows = []
    with open(path) as table:
        for line in table:
            if line.strip():
                rows.append([float(v) for v in line.split()])
    if len(rows) != len(expected_energies):
        raise ValueError("range.txt has %d rows, expected %d (one per /run/beamOn)"
                         % (len(rows), len(expected_energies)))
    for row, energy in zip(rows, expected_energies):
        if len(row) != 7:
            raise ValueError("Expected 7 columns in range.txt, got %r" % row)
        if not math.isclose(row[0] * 1e-6, energy, rel_tol=1e-5):
            raise ValueError("range.txt energy %g eV does not match the macro (%g MeV)"
                             % (row[0], energy))
    return rows


def parse(job):
    rows = read_table(os.path.join(job["path"], "range.txt"), job["ENERGIES"])
    particle = job["PARTICLE"]
    energies = job["ENERGIES"]
    parameters = [{"names": "tracking cut", "values": job["TRACKING_CUT"]}]

    for name, mean, rms in QUANTITIES:
        common = dict(observableName=name,
                      targetName="water",
                      beamParticle=particle,
                      secondaryParticle="None",
                      parameters=parameters,
                      title="%s of %s in liquid water" % (name.capitalize(), particle),
                      xAxisName="E, MeV",
                      yAxisName="%s, nm" % name)
        errors = [row[rms] / math.sqrt(n) for row, n in zip(rows, job["NEVENTS"])]
        yield getJSON(job, "chart",
                      mctool_name="GEANT4",
                      mctool_model=job["PHYSLIST"],
                      beamEnergies=energies,
                      xValues=energies,
                      yValues=[row[mean] for row in rows],
                      yStatErrorsPlus=errors,
                      yStatErrorsMinus=errors,
                      **common)

        if name != "track length":
            continue
        covered = [e for e in energies
                   if e >= icru90.RANGE_MIN_MEV[particle] and icru90.covers(particle, e)]
        if covered:
            # g/cm2 / (g/cm3) = cm = 1e7 nm
            reference = [icru90.csda_range(particle, e) / DENSITY * 1e7 for e in covered]
            yield getJSON(job, "chart",
                          mctool_name="ICRU",
                          mctool_model="ICRU90 CSDA range",
                          beamEnergies=covered,
                          xValues=covered,
                          yValues=reference,
                          **common)
