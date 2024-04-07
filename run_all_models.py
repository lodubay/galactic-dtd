"""
This script runs all multi-zone models for the paper.
"""

from pathlib import Path
import argparse
# Append scripts directory to path to access multizone directory
import sys
import os
sys.path.append(os.path.abspath('./src/scripts'))
import multizone
import paths
import _globals

# Model parameters to be used if others are not specified
_DEFAULT_PARAMS_ = {
    'migration': 'gaussian',
    'evolution': 'insideout',
    'RIa': 'powerlaw',
    'RIa_params': {'slope': -1.1},
    'minimum_delay': _globals.MIN_RIA_DELAY,
    'yields': 'J21',
    'seed': _globals.RANDOM_SEED
}
SFH_LIST = ['insideout', 'lateburst', 'earlyburst', 'twoinfall']
DTD_LIST = ['prompt', 'powerlaw', 'powerlaw', 'exponential', 'exponential',
            'plateau', 'plateau', 'triple']
DTD_PARAMS = [{'peak': 0.05, 'stdev': 0.015, 'timescale': 3.0},
              {'slope': -1.4}, {'slope': -1.1},
              {'timescale': 1.5}, {'timescale': 3.0},
              {'width': 0.3, 'slope': -1.1},
              {'width': 1.0, 'slope': -1.1},
              {'early_rate': 0.05, 'rise_time': 0.5, 'width': 0.5, 'slope': -1.1}]
DTD_NAMES = ['prompt', 'powerlaw_slope14', 'powerlaw_slope11', 
             'exponential_timescale15', 'exponential_timescale30',
             'plateau_width03', 'plateau_width10', 'triple']

def main(seed=_globals.RANDOM_SEED):
    r"""
    Runs the script.
    """
    # Run all combinations of 4 SFHs + 8 DTDs
    count = 1
    total = len(SFH_LIST) * len(DTD_LIST) + 1
    for i, sfh in enumerate(SFH_LIST):
        for j in range(len(DTD_LIST)):
            migration = 'gaussian'
            # Create output directory
            relpath = Path(migration, sfh, DTD_NAMES[j], 'diskmodel')
            fullpath = paths.simulation_outputs / relpath
            if not fullpath.parents[0].exists():
                fullpath.parents[0].mkdir(parents=True)
            # Run counter
            print(f'\nModel {count}/{total}: {str(relpath.parents[0])}')
            params = {
                'migration': migration,
                'evolution': sfh,
                'RIa': DTD_LIST[j],
                'RIa_params': DTD_PARAMS[j],
                'seed': seed,
            }
            model = get_model(str(fullpath), params)
            # model.run([t * model.dt for t in range(round(
            #     _globals.END_TIME / model.dt) + 1)],
            #     overwrite = True, pickle = False)
            count += 1
    # Run fiducial model with analog migration scheme
    relpath = Path('diffusion', 'insideout', 'powerlaw_slope11', 'diskmodel')
    fullpath = paths.simulation_outputs / relpath
    if not fullpath.parents[0].exists():
        fullpath.parents[0].mkdir(parents=True)
    print(f'\nModel {count}/{total}: {str(relpath.parents[0])}')
    model = get_model(str(fullpath), {'migration': 'diffusion', 'seed': seed})
    # model.run([t * model.dt for t in range(round(
    #     _globals.END_TIME / model.dt) + 1)],
    #     overwrite = True, pickle = False)


def get_model(name, params):
    r"""
    Get the milkyway object corresponding to the desired simulation.

    Parameters
    ----------
    name : str
        Path to the simulation output directory, excluding the .vice extension.
    params : dict
        Simulation parameters from Snakemake. The keyword options are:
            
            - 'migration': str, the stellar migration model
            - 'evolution': str, the star formation history model
            - 'RIa': str, the SN Ia delay time distribution model
            - 'RIa_params': dict, kwargs passed to the DTD model
            - 'minimum_delay': float, the minimum SN Ia delay time in Gyr
            - 'yields': str, the nucleosynthetic yield set
            - 'seed': int, seed for random number generator
        
        Any of these parameters left unspecified will be filled by the defaults
        from above.
        
    Returns
    -------
    multizone.src.diskmodel instance
    """
    config = multizone.src.config(
        timestep_size = _globals.DT,
        star_particle_density = _globals.NSTARS,
        zone_width = _globals.ZONE_WIDTH,
        elements = _globals.ELEMENTS
    )
    # Fill in missing params with defaults
    for key, default_value in _DEFAULT_PARAMS_.items():
        if key not in params.keys():
            params[key] = default_value
    kwargs = dict(
        name = name,
        spec = params['evolution'],
        RIa = params['RIa'],
        RIa_kwargs = params['RIa_params'],
        delay = params['minimum_delay'],
        yields = params['yields'],
        seed = params['seed']
    )
    if params['migration'] == 'post-process':
        kwargs['simple'] = True
    else:
        kwargs['migration_mode'] = params['migration']
    return multizone.src.diskmodel.from_config(config, **kwargs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='run_all_models.py',
        description='Run all multi-zone models for the paper.',
    )
    parser.add_argument('-s', '--seed', 
                        metavar='S',
                        type=int,
                        default=_globals.RANDOM_SEED,
                        help='Seed for random number generators.')
    args = parser.parse_args()
    main(**vars(args))
