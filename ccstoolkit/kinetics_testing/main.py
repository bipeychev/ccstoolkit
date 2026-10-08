
## The code has 2 main functions:
	## simulates a kinetic scheme based on user inputted data
	## optimises kinetic parameters (k or A & Ea) form experimental data

from scheme import SCHEME
from simulate import simulate_scheme
from optimise import optimise_parameters


template1 = f"""

[CONDITIONS]
REACTOR TYPE = batch
QSS = false
REACTION TIME = 6000
GRID POINTS = 1000
REACTION TYPE = elementary
TEMPERATURE (K) = 373
PARAMETER TYPE = rate constant

[SPECIES]
RO2R | 4   | stable
RO   | 0   | reactive intermediate
RH   | 8e3 | stable
ROH  | 0   | stable
R    | 0   | reactive intermediate
O2   | 100 | buffered
RO2  | 0   | reactive intermediate
RO2H | 0   | stable

[REACTIONS]
i  | 6.8e-7  | RO2R -> 2 RO
b  | 1.159e9 | RO + RH -> ROH + R
o  | 1e12    | R + O2 -> RO2
a  | 7.15e2  | RO2 + RH -> RO2H + R
dt | 8.55e7  | 2 RO2 -> 2 RO + O2
t  | 1.27e7  | 2RO2 -> RO2R + O2

[EXPERIMENTAL]
DATA TYPE =
DATA = 

"""

template2 = f"""

[CONDITIONS]
REACTOR TYPE = batch
QSS = false
REACTION TIME = 6000
GRID POINTS = 1000
REACTION TYPE = elementary
TEMPERATURE (K) = 373
PARAMETER TYPE = rate constant

[SPECIES]
RO2R | 4   | stable
RO   | 0   | reactive intermediate
RH   | 8e3 | stable
ROH  | 0   | stable
R    | 0   | reactive intermediate
O2   | 100 | buffered
RO2  | 0   | reactive intermediate
RO2H | 0   | stable
Me2CO | 0 | stable
Me | 0 |  reactive intermediate
MeO2 | 0 |  reactive intermediate
MeO2H | 0 | stable

[REACTIONS]
i  | 6.8e-7  | RO2R -> 2 RO
b  | 1.159e9 | RO + RH -> ROH + R
o  | 1e12    | R + O2 -> RO2
a  | 7.15e2  | RO2 + RH -> RO2H + R
dt | 8.55e7  | 2 RO2 -> 2 RO + O2
t  | 1.27e7  | 2RO2 -> RO2R + O2
d  | 6.8e3   | RO -> Me2CO + Me
o2 | 1e12    | Me + O2 -> MeO2
a2 | 6.38e1  | MeO2 + RH -> MeO2H + R

[EXPERIMENTAL]
DATA TYPE =
DATA = 

"""

sim_inputData = template2
#inputData = "input_file_example.txt"


# ----------------- SIMULATE REACTION SCHEME ------------------

#system, results = simulate_scheme(sim_inputData) #,rate_constants) #rate_constants = [1,1,1,1,1,1]

#print(system.scheme._stoichiometric_matrix())
#print(results.raw_data())
#results.print_info()
#results.plot_results()
#results.plot_scheme(filename=None, show=True)
#results.plot_seperate_results()

## info = input (elementary reactions, calculation time, reactor type) + output (stoichiometric matrix, overall reactions,  final concentrations)
#file_name = "test"
#results.export_results(file_name, save_plot=False, save_raw_data=False, save_info=False, save_scheme=False)

# -----------------------------------------------------------------

template = """ 

[CONDITIONS]
REACTOR TYPE = batch
QSS = false
REACTION TIME = {}
GRID POINTS = 100
REACTION TYPE = elementary
TEMPERATURE (K) = {}
PARAMETER TYPE = arrhenius

[REACTIONS]
r1 | 1e6 | 50000 | N2O5 -> NO2 + NO3
r2 | 1 | 1 | NO3 -> NO + O2
r3 | 1 | 1 | NO + NO3 -> 2 NO2 

[SPECIES]
N2O5 | {} | stable
NO2 | 0 | stable
NO3 | 0 | stable
NO | 0 | stable
O2 | 0 | stable

[EXPERIMENTAL]
DATA TYPE = time course
DATA = {}

"""

time2 = 100
time3 = 100
time4 = 100

temp2 = 300
temp3 = 320
temp4 = 340

init2 = 100
init3 = 120
init4 = 150

data2 = {
    "N2O5": {
        "time": [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
        "values": [100.00,98.06,96.08,94.19,  92.31, 90.47, 88.68, 86.93, 85.20, 83.52,81.75]
    }
}

data3 = {
    "N2O5": {
        "time": [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
        "values": [120.00, 112.02, 104.60,97.31,90.55,84.17,78.35,72.86,67.89, 63.15,59.85]
    }
}

data4 = {
    "N2O5": {
        "time": [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
        "values": [150.00, 121.50, 98.60,80.30, 65.40, 53.10,43.10,35.00, 28.40, 23.10,20.15]
    }
}

inputData2 = template.format(time2,temp2,init2,data2)
inputData3 = template.format(time3,temp3,init3,data3)
inputData4 = template.format(time4,temp4,init4,data4)

opt_inputData = [inputData2, inputData3, inputData4]

initA = 1e6
initEa = 50000
initk = 1e3



# ---------------------- OPTIMISE PARAMETERS ----------------------

opt_result = optimise_parameters(opt_inputData, global_search=True) #, optimised_reactions=None, initial_guesses=None, scan=False


# -------------------------------------------------------------------

"""
[REACTIONS]
r3 | -0.86    | NO3 -> NO + O2
r1 | 7.8e-5 | 1e6 | N2O5 <=> NO2 + NO3
r4 | 6e7     | NO3 + NO -> 2NO2
"""





inputData = f"""

[CONDITIONS]
REACTOR TYPE = batch
QSS = false
REACTION TIME = 600
GRID POINTS = 100
REACTION TYPE = net

[REACTIONS]

r1 | k*[N2O5] | 1e-3 | N2O5 -> NO2 + NO3
r2 | k*[NO3] | 1 |  NO3 -> NO  + O2
r3 | k*[NO]*[NO3] | 1 | NO + NO3 -> 2 NO2 

[SPECIES]
N2O5 | 10 | stable
NO2 | 0 | stable
NO3 | 0 | reactive intermediate
NO | 0 |  reactive intermediate
O2 | 0 | stable

"""



inputData = f"""

[CONDITIONS]
REACTOR TYPE = batch
QSS = false
REACTION TIME = 600
GRID POINTS = 100
REACTION TYPE = elementary

[SPECIES]
N2O5 | 10 | stable
NO2 | 0 | stable
NO3 | 0 | reactive intermediate
NO | 0 |  reactive intermediate
O2 | 0 | stable

[REACTIONS]

r1 | 1e-3 | N2O5 -> NO2 + NO3
r2 | 1 |  NO3 -> NO  + O2
r3 | 1 | NO + NO3 -> 2 NO2 

[EXPERIMENTAL]
helooooo

"""




template = """ 

[CONDITIONS]
REACTOR TYPE = batch
QSS = false
REACTION TIME = {}
GRID POINTS = 100
REACTION TYPE = elementary
TEMPERATURE (K) = {}
PARAMETER TYPE = rate constant

[REACTIONS]
r1 | 1e3 | N2O5 -> NO2 + NO3
r2 | 1 | NO3 -> NO + O2
r3 | 1 | NO + NO3 -> 2 NO2 

[SPECIES]
N2O5 | {} | stable
NO2 | 0 | stable
NO3 | 0 | stable
NO | 0 | stable
O2 | 0 | stable

[EXPERIMENTAL]
DATA TYPE = time course
DATA = {}
"""