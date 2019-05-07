def config_parser():
    """Reads in config file and parses it to return important setup configuration parameters"""
    local_dir = None
    ref_spec_dir = None
    ref_fig_size = None
    cal_fig_size = None
    dpi = None
    font_size = None
    imgSizeX = None
    imgSizeY = None

    # with open('C:\\PyCamUV\\config.txt','r') as f:
    with open('config.txt','r') as f:
        for line in f:
            if 'local_dir=' in line:
                local_dir = line.split("'")[1]
            elif 'ref_spec_dir=' in line:
                ref_spec_dir = line.split("'")[1]
            elif 'ref_fig_size=' in line:
                ref_fig_size = line.replace('ref_fig_size=', '').split(' ')[0].split(',')
                ref_fig_size = (float(ref_fig_size[0]), float(ref_fig_size[1]))
            elif 'cal_fig_size=' in line:
                cal_fig_size = line.replace('cal_fig_size=', '').split(' ')[0].split(',')
                cal_fig_size = (float(cal_fig_size[0]), float(cal_fig_size[1]))
            elif 'dpi=' in line:
                dpi = int(line.replace('dpi=', '').split(' ')[0])
            elif 'font_size=' in line:
                font_size = int(line.replace('font_size=', '').split(' ')[0])
            elif 'img_dim=' in line:
                dim = line.replace('img_dim=', '').split(' ')[0].split(',')
                imgSizeX = int(dim[0])
                imgSizeY = int(dim[1])


    return local_dir, ref_spec_dir, ref_fig_size, cal_fig_size, dpi, font_size, imgSizeX, imgSizeY

