from ofe_feh_compare import main

for evolution in ['insideout', 'lateburst']:
    for RIa in ['powerlaw', 'long_delay', 'powerlaw_steep', 'exponential', 'bimodal']:
        output_name = '/'.join(('diffusion', evolution, RIa))
        print('Plotting %s' % output_name)
        main(output_name)