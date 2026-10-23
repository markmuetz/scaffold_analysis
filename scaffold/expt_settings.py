from matplotlib.colors import TABLEAU_COLORS
cycle = list(TABLEAU_COLORS.values())


# Fixed colour per simple experiment, used in every figure (hue = shear, shade = surface
# wind). Checked with a colour-blindness validator; S4W0F's light orange has <3:1 contrast, so
# always give it a legend or direct label.
EXPT_COLOURS = {
    'S0W0Forced': '#4a90d9',
    'S0W5Forced': '#1c5a99',
    'S4W0Forced': '#e08a3c',
    'S4W5Forced': '#a34a12',
}

# Short name, colour, linestyle.
EXPT_DETAILS = {
    'S0W0Forced': ['S0W0F', EXPT_COLOURS['S0W0Forced'], '-'],
    'S4W0Forced': ['S4W0F', EXPT_COLOURS['S4W0Forced'], '-'],
    'S0W5Forced': ['S0W5F', EXPT_COLOURS['S0W5Forced'], '-'],
    'S4W5Forced': ['S4W5F', EXPT_COLOURS['S4W5Forced'], '-'],
    'S0W0rel2_S0W0Forced': ['S0W0_r2_S0W0F', EXPT_COLOURS['S0W0Forced'], '-'],
    'S0W0rel2_S4W0Forced': ['S0W0_r2_S4W0F', EXPT_COLOURS['S0W0Forced'], '--'],
    'S4W0rel2_S0W0Forced': ['S4W0_r2_S0W0F', EXPT_COLOURS['S4W0Forced'], '--'],
    'S4W0rel2_S4W0Forced': ['S4W0_r2_S4W0F', EXPT_COLOURS['S4W0Forced'], '-'],
    'S0W5rel2_S0W5Forced': ['S0W5_r2_S0W5F', EXPT_COLOURS['S0W5Forced'], '-'],
    'S0W5rel2_S4W5Forced': ['S0W5_r2_S4W5F', EXPT_COLOURS['S0W5Forced'], '--'],
    'S4W5rel2_S0W5Forced': ['S4W5_r2_S0W5F', EXPT_COLOURS['S4W5Forced'], '--'],
    'S4W5rel2_S4W5Forced': ['S4W5_r2_S4W5F', EXPT_COLOURS['S4W5Forced'], '-'],
}

for i in range(10):
    EXPT_DETAILS['RWP_C{}'.format(i + 1)] = ['C{}'.format(i + 1), cycle[i], '--']
