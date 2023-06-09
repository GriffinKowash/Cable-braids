# readable version of dscript for reference. Must be run from Discovery.

OUTPUT_DIR = 'C:\\Users\\griffin.kowash\\AppData\\Local\\Temp\\braids_spline_test'
OUTPUT_FILE = 'spline_data.csv'

curve = Selection.GetActive().GetItems[object]()
geometry = curve[0].Shape.RawSubject.Geometry

curve_type = 'spline'  # SplineCurve class is not recognized by API for some reason

if curve_type == 'spline':
    nurbs_data, control_points = curve[0].Shape.RawSubject.Geometry.GetNurbsData()

    degree = nurbs_data.Order - 1
    control_points = [point.Position.ToString()[1:-1] for point in control_points]

    knot_parameters = [knot.Parameter for knot in nurbs_data.Knots]
    knot_multiplicities = [knot.Multiplicity for knot in nurbs_data.Knots]

    knot_vector = []
    for i in range(len(knot_parameters)):
        for _ in range(knot_multiplicities[i]):
            knot_vector.append(knot_parameters[i])

    print(degree)
    print(control_points)
    print(knot_vector)

    file = open(OUTPUT_DIR + '\\' + OUTPUT_FILE, 'w')

    file.write('Degree: ' + str(degree) + '\n')
    file.write('Knot vector: ' + str(knot_vector)[1:-1] + '\n')

    for point in control_points:
        file.write(point + '\n')

    file.close()
    
else:
    print('Type ' + curve_type + ' is not yet supported.')