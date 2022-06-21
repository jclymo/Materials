import pymatgen.core as pg
from z3 import *
import itertools

class z3ConstraintBuilder:
  def __init__(self):
    self.solver = Solver()
    self.variables = {}
  
  def _quantity_variables(self):
    Quantities = self.variables.get('Quantities')
    if Quantities is None:
      Quantities = {f'{str(ion)}': Int(f'n_{str(ion)}') for ion in self.ions}
      self.variables['Quantities'] = Quantities
    return Quantities  
  
  def _element_positive_variables(self):
    Element_Positive = self.variables.get('Element_Positive')
    if Element_Positive is None:
      Element_Positive = {str(element): Bool(f'pos_{str(element)}') for element in self.elements}
      self.variables['Element_Positive'] = Element_Positive
    return Element_Positive

  def _element_negative_variables(self):
    Element_Negative = self.variables.get('Element_Negative')
    if Element_Negative is None:
      Element_Negative = {str(element): Bool(f'neg_{str(element)}') for element in self.elements}
      self.variables['Element_Negative'] = Element_Negative
    return Element_Negative
  
  def _element_count_variables(self):
    Element_Counts = self.variables.get('Element_Counts')
    if Element_Counts is None:
      Element_Counts = {str(element): Int(f'm_{str(element)}') for element in self.elements}
      self.variables['Element_Counts'] = Element_Counts
    return Element_Counts
       
  def non_negativity_constraints(self):
    Quantities = self._quantity_variables() 

    self.solver.add(
      And([X >= 0 for X in Quantities.values()])
    )      

  def charge_constraints(self):        
    Quantities = self._quantity_variables()

    self.solver.add(
      Sum([Quantities[str(ion)]*int(ion.oxi_state) for ion in self.ions]) == 0
    )

  def electronegativity_constraints(self):
    Element_Positive = self._element_positive_variables()
    Element_Negative = self._element_negative_variables()
    Quantities = self._quantity_variables()

    for ion in self.ions:
      if ion.oxi_state > 0:
        self.solver.add(
          Implies(
            Quantities[str(ion)] > 0, 
            Element_Positive[str(ion.element)]
          )
        )
      if ion.oxi_state < 0:
        self.solver.add(
          Implies(
            Quantities[str(ion)] > 0,
            Element_Negative[str(ion.element)]
          )
        )

    for element in self.elements:
      self.solver.add(
        Not(
          And(
            Element_Negative[str(element)], 
            Element_Positive[str(element)]
          )
        )
      )
        
    for element_1 in self.elements:
      for element_2 in self.elements:
        # X is electronegativity (pymatgen)
        if element_1.X < element_2.X:
          self.solver.add(
            Implies(
              Element_Negative[str(element_1)],
              Not(Element_Positive[str(element_2)])
            )
          )

  def element_quantity_constraints(self):
    Quantities = self._quantity_variables()
    Element_Counts = self._element_count_variables()
                
    for element, count in self.element_quantities.items():
      self.solver.add(
        Element_Counts[str(element)] == count
      )

    for element, group in self.element_ions.items():
      self.solver.add(
          Sum([Quantities[str(ion)] for ion in group]) == Element_Counts[str(element)]
      )

  def exclude_solution(self, solution):
    Quantities = self._quantity_variables()
    for multiplier in range(1, self.n_atoms_upper_bound):
      exclusions = []
      for ion_id, quantity in solution.items():
        exclusions.append(Quantities[ion_id] == multiplier*quantity)

      self.solver.add(
        Not(And(*exclusions))
      )


class z3ChargeForCompound(z3ConstraintBuilder):
  def __init__(self, compound, permitted_oxi_states='common'):
    super().__init__()
    self._set_parameters(compound, permitted_oxi_states)
    self._build_constraints()

  def _set_parameters(self, compound, permitted_oxi_states):
    self.n_atoms_upper_bound = int(compound.num_atoms)

    self.elements = compound.elements

    self.element_quantities = dict(compound)

    self.element_ions = {}
    for element in self.elements:
      if permitted_oxi_states == "all":
        self.element_ions[element] = [pg.Species(element.symbol, ch) for ch in element.icsd_oxidation_states or element.oxidation_states]
      elif permitted_oxi_states == "common":
        self.element_ions[element] = [pg.Species(element.symbol, ch) for ch in element.common_oxidation_states]

    self.ions = list(itertools.chain(*self.element_ions.values()))

  def _build_constraints(self):
    self.non_negativity_constraints()
    self.charge_constraints()
    self.electronegativity_constraints()
    self.element_quantity_constraints()
    
  def get_all(self, max_results = 100):
    solutions = []
    Quantities = self._quantity_variables()
    while self.solver.check() == sat and len(solutions) < max_results:
      model = self.solver.model()
      solution = {}
      for ion in self.ions:
        solution[str(ion)] = model[Quantities[str(ion)]].as_long()
      self.exclude_solution(solution)
      solutions.append(solution)
        
    return solutions


def enumerate_and_print_oxidation_states(composition, permitted_oxidation_states=None):
  c = pg.Composition(composition)
  solutions = z3ChargeForCompound(c, permitted_oxidation_states).get_all()
  if len(solutions) == 0:
  	print("No solutions found.")
  for i, solution in enumerate(solutions):
    compact_solution = ', '.join([f'{id}: {num}' for id, num in solution.items() if num > 0])
    print(f"({i}) {compact_solution}")






