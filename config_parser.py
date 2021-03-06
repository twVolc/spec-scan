def config_parser():
    """Reads in config file and parses it to return important setup configuration parameters"""
    # Store all settings in a dictionary, this way I can introduce new settings easily
    config = dict()

    # with open('C:\\PyCamUV\\config.txt','r') as f:
    with open('config.txt', 'r') as f:
        for line in f:
            if 'init_dir=' in line:
                config['init_dir'] = line.split("'")[1]
            elif 'spec_fig_size=' in line:
                spec_fig_size = line.replace('spec_fig_size=', '').split(' ')[0].split(',')
                config['spec_fig_size'] = (float(spec_fig_size[0]), float(spec_fig_size[1]))
            elif 'doas_fig_size=' in line:
                doas_fig_size = line.replace('doas_fig_size=', '').split(' ')[0].split(',')
                config['doas_fig_size'] = (float(doas_fig_size[0]), float(doas_fig_size[1]))
            elif 'scan_fig_size=' in line:
                scan_fig_size = line.replace('scan_fig_size=', '').split(' ')[0].split(',')
                config['scan_fig_size'] = (float(scan_fig_size[0]), float(scan_fig_size[1]))
            elif 'ref_spec_dir=' in line:
                config['ref_spec_dir'] = line.split("'")[1]
            elif 'ref_fig_size=' in line:
                ref_fig_size = line.replace('ref_fig_size=', '').split(' ')[0].split(',')
                config['ref_fig_size'] = (float(ref_fig_size[0]), float(ref_fig_size[1]))
            elif 'cal_fig_size=' in line:
                cal_fig_size = line.replace('cal_fig_size=', '').split(' ')[0].split(',')
                config['cal_fig_size'] = (float(cal_fig_size[0]), float(cal_fig_size[1]))
            elif 'ILS_fig_size=' in line:
                ILS_fig_size = line.replace('ILS_fig_size=', '').split(' ')[0].split(',')
                config['ILS_fig_size'] = (float(ILS_fig_size[0]), float(ILS_fig_size[1]))
            elif 'dpi=' in line:
                config['dpi'] = int(line.replace('dpi=', '').split(' ')[0])
            elif 'font_size_figs=' in line:
                config['mtplt_font_size'] = int(line.replace('font_size_figs=', '').split(' ')[0])
            elif 'arduino_COM' in line:
                config['arduino_COM'] = line.split("'")[1]
            elif 'ILS=' in line:
                config['ILS'] = line.split("'")[1]

            # Reference spectra handling
            elif 'fit_species=' in line:
                species_str = line.split("'")[1]
                config['species'] = species_str.split(',')
            if 'species' in config:
                for species in config['species']:
                    species_id = 'ref_spec_{}'.format(species)
                    if species_id + '=' in line:
                        # Only add this if there is some length of string between the apostrophes in the config file
                        filename = line.split("'")[1]
                        if len(filename) > 0:
                            config[species_id] = filename

        # Perform check that all anticipated reference spectra have been provided. If not, they are deleted from species
        # list and a warning is given
        for species in config['species']:
            species_id = 'ref_spec_{}'.format(species)
            if species_id not in config:
                config['species'].remove(species)
                print('Warning!! {} reference spectra was requested but not provided. It has been removed from processing'.format(species))

    return config

