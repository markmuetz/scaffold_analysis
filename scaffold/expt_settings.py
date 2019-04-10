from matplotlib.colors import TABLEAU_COLORS
cycle = list(TABLEAU_COLORS.values())


# Short name, colour, linestyle.
EXPT_DETAILS = {
    'S0W0Forced': ['S0W0F', cycle[0], '-'],
    'S4W0Forced': ['S4W0F', cycle[1], '-'],
    'S0W5Forced': ['S0W5F', cycle[2], '-'],
    'S4W5Forced': ['S4W5F', cycle[3], '-'],
    'S0W0rel2_S0W0Forced': ['S0W0_r2_S0W0F', cycle[0], '-'],
    'S0W0rel2_S4W0Forced': ['S0W0_r2_S4W0F', cycle[0], '--'],
    'S4W0rel2_S0W0Forced': ['S4W0_r2_S0W0F', cycle[1], '--'],
    'S4W0rel2_S4W0Forced': ['S4W0_r2_S4W0F', cycle[1], '-'],
    'S0W5rel2_S0W5Forced': ['S0W5_r2_S0W5F', cycle[2], '-'],
    'S0W5rel2_S4W5Forced': ['S0W5_r2_S4W5F', cycle[2], '--'],
    'S4W5rel2_S0W5Forced': ['S4W5_r2_S0W5F', cycle[3], '--'],
    'S4W5rel2_S4W5Forced': ['S4W5_r2_S4W5F', cycle[3], '-'],
}

for i in range(10):
    EXPT_DETAILS['RWP_C{}'.format(i + 1)] = ['C{}'.format(i + 1), cycle[i], '-']
