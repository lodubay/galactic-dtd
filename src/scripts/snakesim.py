"""
This script runs a VICE multizone simulation from snakemake parameters.
"""

import multizone
import paths
import _globals

_DEFAULT_PARAMS_ = {
    'migration': 'gaussian',
    'evolution': 'insideout',
    'RIa': 'powerlaw',
    'RIa_params': {'slope': -1.1},
    'minimum_delay': _globals.MIN_RIA_DELAY,
    'yields': 'JW20'
}

def main():
    r"""
    Runs the script.
    """
    fullpath = paths.root / snakemake.output[0] / 'diskmodel'
    if not fullpath.parents[0].exists():
        fullpath.parents[0].mkdir(parents=True)
    model_ = model(str(fullpath), snakemake.params)
    model_.run([_ * model_.dt for _ in range(round(
        _globals.END_TIME / model_.dt) + 1)],
        overwrite = True, pickle = False)


def model(name, params):
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
    for key, value in _DEFAULT_PARAMS_.items():
        if key not in params.keys():
            params[key] = value
    kwargs = dict(
        name = name,
        spec = params.evolution,
        RIa = params.RIa,
        RIa_kwargs = params.RIa_params,
        delay = params.minimum_delay,
        yields = params.yields
    )
    if params.migration == "post-process":
        kwargs["simple"] = True
    else:
        kwargs["migration_mode"] = params.migration
    return multizone.src.diskmodel.from_config(config, **kwargs)


if __name__ == '__main__':
    main()
