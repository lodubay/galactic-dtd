r"""
This file contains a collection of colorblind-safe color palettes by Paul Tol.
More information can be found at the website https://personal.sron.nl/~pault/
"""

from matplotlib.colors import ListedColormap

###############################################################################
#                                                                             #
#       QUALITATIVE COLORMAPS                                                 #
#                                                                             #
###############################################################################

# Bright: for lines and labels
bright = ListedColormap(
    ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377',
     '#BBBBBB'],
    name="paultol_bright")
# High-contrast: for monochrome vision and printouts
highcontrast = ListedColormap(
    ['#004488', '#DDAA33', '#BB5566'],
    name='paultol_highcontrast')
# Vibrant: for TensorBoard
vibrant = ListedColormap(
    ['#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988',
     '#BBBBBB'],
    name='paultol_vibrant')
# Muted: for more colors
muted = ListedColormap(
    ['#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255',
     '#44AA99', '#999933', '#AA4499', '#DDDDDD'],
    name="paultol_muted")
# Medium-contrast: for color pairs (1 with 4, 2 with 5, 3 with 6)
mediumcontrast = ListedColormap(
    ['#6699CC', '#004488', '#EECC66', '#994455', '#997700', '#EE99AA'],
    name="paultol_mediumcontrast")
# Pale: for the background of black text
pale = ListedColormap(
    ['#BBCCEE', '#CCEEFF', '#CCDDAA', '#EEEEBB', '#FFCCCC', '#DDDDDD'],
    name="paultol_pale")
# Dark: for the background of white text or for text on a white background
dark = ListedColormap(
    ['#222255', '#225555', '#225522', '#666633', '#663333', '#555555'],
    name="paultol_dark")
# Light: for filling labeled cells
light = ListedColormap(
    ['#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF', '#44BB99',
     '#BBCC33', '#AAAA00', '#DDDDDD'],
    name="paultol_light")
