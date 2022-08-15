r"""
This script generates plots for both the single degenerate and double
degenerate cases from Grettio 2005. Reference this script in the `\script{}`
command when plotting both figures side-by-side in the LaTeX document.
"""

from onezone_greggio05_single import main as main_single
from onezone_greggio05_double import main as main_double

main_single()
main_double()