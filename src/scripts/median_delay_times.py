"""
This script outputs the median delay times of each DTD form.
"""

import numpy as np
from delay_time_distributions import styles
import paths
import _globals

def main():
    distributions = [styles.prompt, styles.plaw_steep, styles.plaw,
                     styles.exp, styles.exp_long, styles.plateau,
                     styles.plateau_long, styles.triple]
    names = [dtd['label'] for dtd in distributions]
    dt = 0.001
    times = np.arange(_globals.MIN_RIA_DELAY, _globals.END_TIME + dt, dt)
    medians = []
    for dtd in distributions:
        func = dtd['func']
        yvals = np.array([func(t) for t in times])
        cdf = np.cumsum(yvals / np.sum(yvals))
        med_idx = np.where(cdf >= 0.5)[0][0]
        medians.append(times[med_idx])
    # Write output file
    with open(paths.output / 'median_delay_times.txt', 'w') as f:
        f.writelines('%s\t%s\n' % (n, round(m, 3)) for n, m in zip(names, medians))


if __name__ == '__main__':
    main()
