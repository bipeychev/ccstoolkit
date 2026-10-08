
## Contains the simulate function used to run the code

import time
from scheme import SCHEME
from reactors import BATCH_REACTOR, PLUG_FLOW_REACTOR, CSTR_REACTOR
from qss_reactors import BATCH_REACTOR_QSS, PLUG_FLOW_REACTOR_QSS, CSTR_REACTOR_QSS
from results import RESULTS

def simulate_scheme(inputData=None, rate_constants=None):
	overall_start_time = time.perf_counter()

	scheme = SCHEME(inputData)

	if scheme._N_grid is not None:
		N_grid = scheme._N_grid
	else:
		N_grid = 1000

	if rate_constants is not None:
		scheme._set_rate_constants(rate_constants)

	reactor = scheme._reactor_type
	qss = scheme._qss_true
	reaction_type = scheme._reaction_type

	REACTORS = {
		"batch": BATCH_REACTOR,
		"pfr": PLUG_FLOW_REACTOR,
		"cstr": CSTR_REACTOR,
	}

	QSS_REACTORS = {
		"batch": BATCH_REACTOR_QSS,
		"pfr": PLUG_FLOW_REACTOR_QSS,
		"cstr": CSTR_REACTOR_QSS,
	}

	if qss:
		reactor_class = QSS_REACTORS[reactor]
	else:
		reactor_class = REACTORS[reactor]
	
	system = reactor_class(scheme)
	
	if reactor == "cstr":
		system.simulate(system.scheme._reaction_time, system.scheme._residence_time, N_grid=N_grid)
	elif reactor == "pfr":
		system.simulate(system.scheme._residence_time, N_grid=N_grid)
	elif reactor == "batch":
		system.simulate(system.scheme._reaction_time, N_grid=N_grid)

	overall_end_time = time.perf_counter()
	system.overall_runtime = overall_end_time - overall_start_time

	results = RESULTS(system, reactor, qss, reaction_type)
	
	return system, results

