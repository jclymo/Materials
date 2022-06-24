import pymatgen.core as pg
import ChargesSolver
import csv
from pathlib import Path
import json

# data file from pymatgen. 
# required here since the interface doesn't seem to have a way 
# of retrieving the list of shanon radii
with open(str(Path(__file__).absolute().parent / "periodic_table.json")) as f:
    _pt_data = json.load(f)

def shannon_oxi_states(element):
  el_radii = _pt_data[element.symbol]['Shannon radii']
  return [int(ch) for ch in el_radii.keys()]

def ions_for_elements(elements, permitted_oxi_states):
  if permitted_oxi_states == "shannon":
    return {element: 
      [
        pg.Species(element.symbol, ch) for ch in shannon_oxi_states(element)
      ] 
      for element in elements
    }
  elif permitted_oxi_states == "common":
    return {element: 
      [
        pg.Species(element.symbol, ch) for ch in element.common_oxidation_states
      ] 
      for element in elements
    }
  elif permitted_oxi_states == "all":
    return {element: 
      [
        pg.Species(element.symbol, ch) for ch in element.icsd_oxidation_states or element.oxidation_states
      ] 
      for element in elements
    }
  raise Exception("Value {permitted_oxi_states} for permitted_oxi_states not recognised.")

def get_oxid_state_guesses(composition, permitted_oxi_states="shannon"):
  c = pg.Composition(composition)
  element_ions = ions_for_elements(c.elements, permitted_oxi_states)
  return ChargesSolver.z3ChargeForCompound(dict(c), element_ions).get_all()
  








