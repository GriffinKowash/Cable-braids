import sys, getopt, json, time, traceback
from braid import Braid

try:
    # Command line parameters
    arg_str = 'm:g:f:o:r:'
    arg_list = ['mode =', 'geo =', 'file =', 'outdir =', 'resolution =']
    options, args = getopt.getopt(sys.argv[1:], arg_str, arg_list)
    
    for opt, arg in options:
        if opt in ['-m', '--mode ']:
            print('mode: |||' + arg + '|||')
            if arg == 'spline':
                mode = arg
            elif arg == 'line':
                print('Line mode is not currently supported via the Discovery API.')
            elif arg == 'points':
                print('Point mode is not currently supported via the Discovery API.')
            else:
                print(f'Mode {arg} not recognized.')  
        elif opt in ['-g', '--geo ']:
            geo_path = arg
        elif opt in ['-f', '--file ']:
            file_path = arg
        elif opt in ['-o', '--outdir ']:
            output_dir = arg
        elif opt in ['-r', '--resolution ']:
            resolution = int(arg)
        else:
            print(f'Unrecognized option {opt} to discobraid.py')
   
    with open(geo_path + '\\geo.json', 'r') as geo_file:
        geo_dict = json.load(geo_file)
    
    braid = Braid()
    braid.set_geometry(geo_dict)
    
    if mode == 'spline':
        with open(file_path + '\\spline.json', 'r') as spline_file:
            spline_dict = json.load(spline_file)
        print('resolution: ', resolution, type(resolution))
        print('spline: ', spline_dict)

        braid.set_path_from_spline(spline_dict, resolution, True)
        
    
    braid.construct(verbose=False)
    braid.save(output_dir)
        
except:
    traceback.print_exc()
    time.sleep(100)