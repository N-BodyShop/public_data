from tangos.properties.pynbody import PynbodyPropertyCalculation
from tangos.properties import PropertyCalculation
import pynbody
import numpy as np


class Outflows(PynbodyPropertyCalculation):
    names = ""
    def calculate(self, particle_data, existing_properties):
       pass
   def requires_property(self):
       return ["shrink_center", "max_radius", "vel_center"]
