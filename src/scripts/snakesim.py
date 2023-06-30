"""
This script runs a VICE multizone simulation from snakemake parameters.
"""

import multizone
import paths
import _globals

def main():
    r"""
    Runs the script.
    """
    fullpath = paths.root / snakemake.output / 'diskmodel'
    model_ = model(str(fullpath), snakemake.params)
    model_.run([_ * model_.dt for _ in range(round(
        _globals.END_TIME / model_.dt) + 1)],
        overwrite = True, pickle = False)


def model(name, params):
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
