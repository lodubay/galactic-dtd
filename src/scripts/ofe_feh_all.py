from ofe_feh_compare import main

for evolution in ['insideout', 'lateburst', 'conroy22']:
    for RIa in ['powerlaw', 'powerlaw_delayed', 'powerlaw_steep',
                'powerlaw_steep_delayed', 'powerlaw_broken', 'exponential',
                'exponential_long', 'exponential_delayed', 'bimodal']:
        output_name = '/'.join(('diffusion', evolution, RIa))
        print('Plotting %s' % output_name)
        main(output_name)