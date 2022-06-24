import pymatgen.core as pg
from z3 import *
import itertools

class z3ConstraintBuilder:
  def __init__(self):
    self.solver = Solver()
    self.variables = {}
  
  def _variables(self, var_name, var_type, ids):
    local_variables = self.variables.get(var_name)
    if local_variables is None:
      local_variables = {str(idx): var_type(f'{var_name}_{str(idx)}') for idx in ids}
      self.variables[var_name] = local_variables
    return local_variables

  def _species_quantity_variables(self):
    return self._variables('Quantities', Int, self.ions)
  
  def _element_positive_variables(self):
    return self._variables('Element_Positive', Bool, self.elements)

  def _element_negative_variables(self):
    return self._variables('Element_Negative', Bool, self.elements)
  
  def _element_count_variables(self):
    return self._variables('Element_Count', Int, self.elements)    
       
  def non_negativity_constraints(self):
    Quantities = self._species_quantity_variables() 

    self.solver.add(
      And([X >= 0 for X in Quantities.values()])
    )      

  def charge_constraints(self):        
    Quantities = self._species_quantity_variables()

    self.solver.add(
      Sum([Quantities[str(ion)]*int(ion.oxi_state) for ion in self.ions]) == 0
    )

  def electronegativity_constraints(self):
    Element_Positive = self._element_positive_variables()
    Element_Negative = self._element_negative_variables()
    Quantities = self._species_quantity_variables()

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
        # X is electronegativity
        if element_1.X < element_2.X:
          self.solver.add(
            Implies(
              Element_Negative[str(element_1)],
              Not(Element_Positive[str(element_2)])
            )
          )

  def element_quantity_constraints(self):
    Quantities = self._species_quantity_variables()
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
    Quantities = self._species_quantity_variables()
    for multiplier in range(1, self.n_atoms_upper_bound):
      exclusions = []
      for ion_id, quantity in solution.items():
        exclusions.append(Quantities[ion_id] == multiplier*quantity)

      self.solver.add(
        Not(And(*exclusions))
      )


class z3ChargeForCompound(z3ConstraintBuilder):
  def __init__(self, element_quantities, element_ions):
    super().__init__()
    
    if any([n != int(n) for n in element_quantities.values()]):
      raise Exception("Elements must have integer quantities.")

    self.element_quantities = {element: int(n) for element, n in element_quantities.items()}
    self.elements = element_quantities.keys()
    self.element_ions = element_ions
    self.ions = list(itertools.chain(*self.element_ions.values()))
    self.n_atoms_upper_bound = sum(self.element_quantities.values())

    self._build_constraints()

  def _build_constraints(self):
    self.non_negativity_constraints()
    self.charge_constraints()
    self.electronegativity_constraints()
    self.element_quantity_constraints()
    
  def get_all(self, max_results=100):
    solutions = []
    Quantities = self._species_quantity_variables()
    while self.solver.check() == sat and len(solutions) < max_results:
      model = self.solver.model()
      solution = {}
      for ion in self.ions:
        solution[str(ion)] = model[Quantities[str(ion)]].as_long()
      self.exclude_solution(solution)
      solutions.append(solution)
        
    return solutions



