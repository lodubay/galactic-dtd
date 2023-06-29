"""
This script runs all necessary multi-zone simulations for the paper.
"""

import simulations
import paths

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

def main():
    dtd = DTD_LIST[0]
    model('insideout', dtd['model'], dtd_kwargs=dtd['params'], dtd_name=dtd['name'])


def model(sfh, dtd, dtd_kwargs={}, dtd_name='', migration='gaussian'):
    r"""
    Get the milkyway object corresponding to the desired simulation.

    Parameters
    ----------
    args : argparse.Namespace
        The command line arguments parsed via argparse.
    """
    config = simulations.config(
        timestep_size = simulations._globals.DT,
        star_particle_density = simulations._globals.NSTARS,
        zone_width = simulations._globals.ZONE_WIDTH,
        elements = simulations._globals.ELEMENTS
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
        delay = simulations._globals.MIN_RIA_DELAY,
        yields = 'JW20'
    )
    if migration == "post-process":
        kwargs["simple"] = True
    else:
        kwargs["migration_mode"] = migration
    return simulations.diskmodel.from_config(config, **kwargs)


if __name__ == '__main__':
    main()
