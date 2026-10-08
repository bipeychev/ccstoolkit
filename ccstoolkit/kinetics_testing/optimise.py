
## contains the optimise_parameters function used for general optimisation

import numpy as np
import matplotlib.pyplot as plt
from scheme import SCHEME
from simulate import simulate_scheme
from arrhenius_fitting import ARRHENIUS_FITTING
from rc_fitting import RC_FITTING

def optimise_parameters(inputData=None, optimised_reactions=None, initial_guesses=None, scan=False, global_search=False):

	if inputData is None:
		raise ValueError(f"Input data must be provided")
	if isinstance(inputData, str):
		inputData = [inputData]
	scheme = SCHEME(inputData[0])

	# optimises for rate constants or arrhenius parameters
	if scheme._parameter_type == "rate constant":
		fitting = RC_FITTING(inputData=inputData, optimised_reactions=optimised_reactions, initial_guesses=initial_guesses)
		result = fitting.optimise_rate_constants(scan=scan)

	elif scheme._parameter_type == "arrhenius":
		fitting = ARRHENIUS_FITTING(inputData=inputData, optimised_reactions=optimised_reactions, initial_guesses=initial_guesses)
		result = fitting.optimise_arrhenius(scan=scan, global_search=global_search)

	return result