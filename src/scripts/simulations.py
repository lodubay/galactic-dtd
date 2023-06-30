"""
This script runs all necessary multi-zone simulations for the paper.
"""

import multizone
import paths
import _globals

# SN Ia delay time distribution models
DTD_LIST = [
    {'name': 'powerlaw_slope11',
     'model': 'powerlaw',
     'params': {'slope': -1.1}},
    {'name': 'powerlaw_slope14',
     'model': 'powerlaw',
     'params': {'slope': -1.4}},
    {'name': 'exponential_timescale15',
     'model': 'exponential',
     'params': {'timescale': 1.5}},
    {'name': 'exponential_timescale30',
     'model': 'exponential',
     'params': {'timescale': 3.}},
    {'name': 'plateau_width300_slope11',
     'model': 'plateau',
     'params': {'width': 0.3, 'slope': -1.1}},
    {'name': 'plateau_width1000_slope11',
     'model': 'plateau',
     'params': {'width': 0.3, 'slope': -1.1}},
    {'name': 'prompt_peak050_stdev015_timescale30',
     'model': 'prompt',
     'params': {'peak': 0.05, 'stdev': 0.015, 'timescale': 3.}},
    {'name': 'triple_rise500_width500_slope11',
     'model': 'triple',
     'params': {'rise_time': 0.5, 'width': 0.5, 'slope': -1.1}},
]
# Star formation history models
SFH_LIST = ['insideout', 'lateburst', 'earlyburst', 'twoinfall']

def main():
    for sfh in SFH_LIST:
        for dtd in DTD_LIST:
            model_ = model('insideout', dtd['model'], dtd_kwargs=dtd['params'], 
                           dtd_name=dtd['name'])
            model_.run([_ * model_.dt for _ in range(round(
                _globals.END_TIME / model_.dt) + 1)],
                overwrite = False, pickle = False)


def model(sfh, dtd, dtd_kwargs={}, dtd_name='', migration='gaussian'):
    r"""
    Get the milkyway object corresponding to the desired simulation.

    Parameters
    ----------
    args : argparse.Namespace
        The command line arguments parsed via argparse.
    """
    config = multizone.src.config(
        timestep_size = _globals.DT,
        star_particle_density = _globals.NSTARS,
        zone_width = _globals.ZONE_WIDTH,
        elements = _globals.ELEMENTS
    )
    # Path to output directory
    if dtd_name == '':
        dtd_name = dtd
    fullpath = paths.simulation_outputs / migration / sfh / dtd_name
    if not fullpath.parents[0].exists():
        fullpath.parents[0].mkdir(parents=True)
    kwargs = dict(
        name = str(fullpath),
        spec = sfh,
        RIa = dtd,
        RIa_kwargs = dtd_kwargs,
        delay = _globals.MIN_RIA_DELAY,
        yields = 'JW20'
    )
    if migration == "post-process":
        kwargs["simple"] = True
    else:
        kwargs["migration_mode"] = migration
    return multizone.src.diskmodel.from_config(config, **kwargs)


if __name__ == '__main__':
    main()
