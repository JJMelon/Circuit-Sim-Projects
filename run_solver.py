from scripts.Solve import solve

# path to the grid network RAW file
casename = 'testcases/ACTIVSg500_prior_solution.RAW'

# the settings for the solver
settings = {
    "Tolerance": 1E-05,
    "Max Iters": 1000,
    "Limiting":  False,
    "Flat Start": False
}

# run the solver
solve(casename, settings)
