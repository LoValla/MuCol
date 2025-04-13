#####################################
#
# simple script to create lcio files with single particle
# events - modify as needed
# @author F.Gaede, DESY
# @date 1/07/2014
#
# initialize environment:
#  export PYTHONPATH=${LCIO}/src/python:${ROOTSYS}/lib
#
#####################################
import numpy as np
import random
from array import array
from g4units import deg, s

# --- LCIO dependencies ---
from pyLCIO import EVENT, IMPL, IOIMPL
# from pyLCIO import UTIL, IO

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--numberOfEvents", type = int, default = 10000)
#parser.add_argument("--gunEnergy", type = float)
parser.add_argument("--output", type = str, default="tau_guns.slcio")
args = parser.parse_args()

# ---- number of events per momentum bin -----
nevt = args.numberOfEvents

outfile = args.output

# --------------------------------------------

wrt = IOIMPL.LCFactory.getInstance().createLCWriter()

wrt.open(outfile, EVENT.LCIO.WRITE_NEW)

random.seed()


# ========== particle properties ===================

# momenta = [ 1. , 3., 5., 10., 15., 25., 50., 100. ]
# momenta = [ 5. ]

# generate flat in E
#E = args.gunEnergy

# for i in range(10):
#    pT.append(random.random()*300.+20.)

genstat = 1
pdg = -15
mass = 1.77686 # tau mass
charge = +1.0

# decay time in seconds
lifetime = 2.9e-13 * s

# bounds on theta and conversion to radians
theta_min = 10.0 * deg
theta_max = 170.0 * deg

# =================================================


for j in range(0, nevt):
    col = IMPL.LCCollectionVec(EVENT.LCIO.MCPARTICLE)
    evt = IMPL.LCEventImpl()

    evt.setEventNumber(j)

    evt.addCollection(col, "MCParticle")

    phi = random.random() * np.pi * 2.0  # flat in phi
    # theta = np.acos(2*random.random()-1) #flat in eta
    theta = random.uniform(theta_min, theta_max)
    
    E = random.betavariate(2., 4.75) * 1000. + mass

    p = np.sqrt(E**2 - mass**2)

    px = p * np.sin(theta) * np.cos(phi)
    py = p * np.sin(theta) * np.sin(phi)
    pz = p * np.cos(theta)

    momentum = array("f", [px, py, pz])

    gamma = E / mass

    beta = p / E

    decayTime = np.random.exponential(lifetime)

    # c in mm/ns
    decayLen = gamma * beta * 299.8 * decayTime

    epx = decayLen * np.cos(phi) * np.sin(theta)
    epy = decayLen * np.sin(phi) * np.sin(theta)
    epz = decayLen * np.cos(theta)

    endpoint = array("d", [epx, epy, epz])

    # --------------- create MCParticle -------------------

    mcp = IMPL.MCParticleImpl()

    mcp.setGeneratorStatus(genstat)
    mcp.setMass(mass)
    mcp.setPDG(pdg)
    mcp.setMomentum(momentum)
    mcp.setCharge(charge)
    # set endpoint (remember G4 has a 0.7mm threshold cut for secondary prts production)
    mcp.setEndpoint(endpoint)

    # -------------------------------------------------------

    col.addElement(mcp)

    wrt.writeEvent(evt)


wrt.close()


#
#  longer format: - use ".hepevt"
#

#
#    int ISTHEP;   // status code
#    int IDHEP;    // PDG code
#    int JMOHEP1;  // first mother
#    int JMOHEP2;  // last mother
#    int JDAHEP1;  // first daughter
#    int JDAHEP2;  // last daughter
#    double PHEP1; // px in GeV/c
#    double PHEP2; // py in GeV/c
#    double PHEP3; // pz in GeV/c
#    double PHEP4; // energy in GeV
#    double PHEP5; // mass in GeV/c**2
#    double VHEP1; // x vertex position in mm
#    double VHEP2; // y vertex position in mm
#    double VHEP3; // z vertex position in mm
#    double VHEP4; // production time in mm/c
#
#    inputFile >> ISTHEP >> IDHEP
#    >> JMOHEP1 >> JMOHEP2
#    >> JDAHEP1 >> JDAHEP2
#    >> PHEP1 >> PHEP2 >> PHEP3
#    >> PHEP4 >> PHEP5
#    >> VHEP1 >> VHEP2 >> VHEP3
#    >> VHEP4;
