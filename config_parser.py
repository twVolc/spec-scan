def config_parser():
    """Reads in config file and parses it to return important setup configuration parameters"""
    # Store all settings in a dictionary, this way I can introduce new settings easily
    config = dict()

    # with open('C:\\PyCamUV\\config.txt','r') as f:
    with open('config.txt','r') as f:
        for line in f:
            if 'local_dir=' in line:
                config['local_dir'] = line.split("'")[1]
            elif 'ref_spec_dir=' in line:
                config['ref_spec_dir'] = line.split("'")[1]
            elif 'ref_fig_size=' in line:
                ref_fig_size = line.replace('ref_fig_size=', '').split(' ')[0].split(',')
                config['ref_fig_size'] = (float(ref_fig_size[0]), float(ref_fig_size[1]))
            elif 'cal_fig_size=' in line:
                cal_fig_size = line.replace('cal_fig_size=', '').split(' ')[0].split(',')
                config['cal_fig_size'] = (float(cal_fig_size[0]), float(cal_fig_size[1]))
            elif 'dpi=' in line:
                config['dpi'] = int(line.replace('dpi=', '').split(' ')[0])
            elif 'font_size_figs=' in line:
                config['mtplt_font_size'] = int(line.replace('font_size_figs=', '').split(' ')[0])
            elif 'ref_spec_SO2=' in line:
                config['ref_spec_SO2'] = line.split("'")[1]
            elif 'ILS=' in line:
                config['ILS'] = line.split("'")[1]
            elif 'arduino_COM' in line:
                config['arduino_COM'] = line.split("'")[1]

    return config

