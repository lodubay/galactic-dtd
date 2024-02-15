"""
This script generates two plots to be included in a single figure: the effect
of the minimum SN Ia delay time and the different DTD models on the one-zone
model outputs.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import vice
import paths
from multizone.src.yields import J21
from multizone.src import models
from _globals import END_TIME, ONEZONE_DEFAULTS, TWO_COLUMN_WIDTH
from track_and_mdf import setup_axes, plot_vice_onezone
from delay_time_distributions import styles

def main(overwrite=False):
    plt.style.use(paths.styles / 'paper.mplstyle')
    fig = plt.figure(figsize=(TWO_COLUMN_WIDTH, 0.54*TWO_COLUMN_WIDTH))
    gs = fig.add_gridspec(7, 15, wspace=0.)
    # subfigs = fig.subfigures(1, 3, wspace=0.)
    subfigs = [fig.add_subfigure(gs[:,i:i+w]) for i, w in zip((0, 8), (8, 7))]
    minimum_delay(subfigs[0])
    dtds(subfigs[1])
    plt.subplots_adjust(bottom=0.09, top=0.98, left=0.11, right=0.98, wspace=0.5)
    fig.savefig(paths.figures / 'onezone_twopanel.pdf', dpi=300)
    plt.close()
    

def minimum_delay(subfig, delays=[0.15, 0.04], line_styles=['--', '-']):
    output_dir = paths.data / 'onezone' / 'minimum_delay'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    axs = setup_axes(subfig)

    dt = ONEZONE_DEFAULTS['dt']
    simtime = np.arange(0, END_TIME + dt, dt)
    params = ONEZONE_DEFAULTS.copy()

    prev_eta = 0.
    for delay, ls in zip(delays, line_styles):
        params['delay'] = delay
        distributions = [styles.exp, styles.plateau, styles.plaw]
        for i, dtd in enumerate(distributions):
            if delay == delays[1]:
                label = dtd['label']
            else:
                label = None
            # Modify mass-loading factor for clarity
            if i == 0:
                params['eta'] = 1.
            else:
                params['eta'] = ONEZONE_DEFAULTS['eta']
            name = dtd['func'].name + '_delay{:03d}'.format(int(delay * 1000))
            sz = vice.singlezone(name=str(output_dir / name),
                                 RIa=dtd['func'],
                                 func=models.insideout(8, dt=dt), 
                                 mode='sfr',
                                 **params)
            sz.run(simtime, overwrite=True)
            plot_vice_onezone(str(output_dir / name), 
                              fig=subfig, axs=axs,
                              label=label, 
                              color=dtd['color'],
                              linestyle=ls,
                              zorder=i,
                              marker_labels=(i==0 and delay==delays[0])
                              )
            # label with eta value
            if ls == '-' and prev_eta != params['eta']:
                data = axs[0].lines[-1].get_xydata()
                textxy = data[-1]
                if i == 1:
                    label = r'$\eta=%s$' % params['eta']
                else:
                    label = str(params['eta'])
                axs[0].text(textxy[0]+0.1, textxy[1]-0.02, label, 
                            va='top', ha='right')
            prev_eta = params['eta']


    axs[0].set_xlim(right=0.7)
    # Re-scale marginal axis limits
    axs[1].set_ylim(bottom=0)
    axs[2].set_xlim(left=0)
    # Legend
    handles, labels = axs[0].get_legend_handles_labels()
    handles += [Line2D([], [], color='k', ls=ls, lw=1) for ls in line_styles]
    labels += [f'{int(1000*delay)} Myr minimum delay' for delay in delays]
    axs[0].legend(handles, labels, frameon=False, loc='lower left',
                  handlelength=1.2)
    return axs


def dtds(subfig):
    output_dir = paths.data / 'onezone' / 'dtd'
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    
    axs = setup_axes(subfig, ylabel=False)
    
    dt = ONEZONE_DEFAULTS['dt']
    simtime = np.arange(0, END_TIME + dt, dt)
    
    distributions = [styles.triple, styles.plateau_long, styles.exp, 
                     styles.plaw, styles.prompt]
    
    for i, dtd in enumerate(distributions):
        sz = vice.singlezone(name=str(output_dir / dtd['func'].name),
                             RIa=dtd['func'],
                             func=models.insideout(8, dt=dt), 
                             mode='sfr',
                             **ONEZONE_DEFAULTS)
        sz.run(simtime, overwrite=True)
        plot_vice_onezone(str(output_dir / dtd['func'].name), 
                          fig=subfig, axs=axs,
                          label=dtd['label'], 
                          color=dtd['color'],
                          linestyle=dtd['line'],
                          marker_labels=(i==0),
                          )
    
    # Re-scale marginal axis limits
    axs[1].set_ylim(bottom=0)
    axs[2].set_xlim(left=0)
    axs[0].legend(frameon=False, loc='lower left', handlelength=1.8)
    return axs


if __name__ == '__main__':
    main()
