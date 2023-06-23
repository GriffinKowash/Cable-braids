# Python Script, API Version = V23
import os
import csv
import time
import string

SOURCE_DIR = 'C:\\Users\\griffin.kowash\\AppData\\Local\\Temp\\braids\\braid_data'
SCRIPT_DIR = 'C:\\Users\\griffin.kowash\\Documents\\Projects\\Cable braids'
SCRIPT_NAME = 'GenerateBraid.py'
SCALE = 0.001  # convert mm to m

load_external = False


####### GEOMETRY PARAMETERS #######

if not load_external:
    s = 1.78 / 2       # shield radius (mm)
    c = 8       # number of carriers (count)
    n = 4       # "ends"/wires per carrier (count)
    d = 0.131   # wire diameter (mm)
    p = 100     # "picks"/crossing frequency (1/mm)
    alpha = 72  # pitch angle (degrees)
    z_max = 145        # endpoint of z axis (mm)
    resolution = 400  # number of points per curve (count)

    plotting = False           # whether to visualize braid in matplotlib
    plot_mode = 'line'     # plot using lines ('line') or surfaces ('surf')

    command = 'python {}\\{} -s %f -c {} -n {} -d {} -a {} -z {} -r {} -o {} -m {} --save'.format(SCRIPT_DIR, SCRIPT_NAME, s, c, n, d, alpha, z_max, resolution, SOURCE_DIR, plot_mode)
    if plotting:
        command = command + ' --plot'
        
    os.system(command)
    
####### END GEOMETRY PARAMTERS #######


# Load config info
config_file = open(SOURCE_DIR + '\\config.txt', 'r')
config_data = csv.reader(config_file, delimiter=',')

for row in config_data:
    print(row)
    label, value = row
    
    if label == 'carriers':
        carriers = int(value)
    elif label == 'wires per carrier':
        wires = int(value)
        
config_file.close()

# Sketch braid curves
Sketch3D.Set3DSketchMode(True)
start_time = time.time()

for chirality in ['cw', 'ccw']:
    print('_______' + chirality + ' direction_______')
    for i in range(carriers // 2):
        print('    CARRIER ' + str(i))
        for j in range(wires):
            print('        Wire ' + str(j))
            data_file = open(SOURCE_DIR + '\\data_' + chirality + '_c' + str(i) + '_w' + str(j) +  '.csv', 'r')
            data = csv.reader(data_file, delimiter=',')
            
            points = []
            for row in data:
                x, y, z = row
                points.append(Point.Create(SCALE * float(x), SCALE * float(y), SCALE * float(z)))
                
            data_file.close()
            SketchLine.CreateChain(points)

time_elapsed = time.time() - start_time
Sketch3D.Set3DSketchMode(False)
#mode = InteractionMode.Solid
#result = ViewHelper.SetViewMode(mode, Info6)

print('Braid import complete.')
print('Elapsed time: ' +str(int(time_elapsed / 60)) + ' minutes and ' + str(round(time_elapsed % 60)) + ' seconds.')