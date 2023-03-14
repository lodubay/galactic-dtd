from ofe_feh_compare import main

for evolution in ['insideout', 'lateburst', 'conroy22', 'twoinfall', 'conroy22_JW20yields']:
    for RIa in ['powerlaw_slope11', 
                'powerlaw_slope14', 
                'exponential_timescale15',
                'exponential_timescale30',
                'plateau_width300_slope11',
                'plateau_width1000_slope11',
                'prompt_peak050_stdev015_timescale30',
                'triple_delay040']:
        output_name = '/'.join(('diffusion', evolution, RIa))
        print('Plotting %s' % output_name)
        main(output_name)