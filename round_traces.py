"""
This script modifies an eagle board file by rounding corners of traces.



Process for creating rounded traces:




"""

import os
import pdb
import math
from xml.etree import ElementTree

from point2d import Point2D

# filename to modify
filename = r'C:\Users\tdkostk\Documents\eagle\projects\round_traces\round_traces_test.brd'
filename = r'C:\Users\tdkostk\Documents\eagle\projects\micro_ohmmeter\micro_ohmmeter_rev5.brd'
filename = r'C:\Users\tdkostk\Documents\eagle\projects\kct-tester\sandia-cable-tester-rev5.brd'


class StraightWire:
    """Holds information about a straight trace."""

    def __init__(self, wire, signal=''):
        assert 'curve' not in wire.attrib
        self.start_point = Point2D(float(wire.attrib['x1']),
                                   float(wire.attrib['y1']))
        self.end_point = Point2D(float(wire.attrib['x2']),
                                 float(wire.attrib['y2']))
        self.length = (self.start_point - self.end_point).norm()
        self.signal = signal
        self.width = wire.attrib['width']
        self.layer = wire.attrib['layer']


# read in wires
if False:
    tree = ElementTree.parse(filename)
    root = tree.getroot()

    # hold wires for each signal
    wire = dict()

    for signal in root.iter('signal'):
        print(signal.tag, signal.attrib)
        # read in all wires
        wires = []
        # for each point, count the number of wires that connect it
        wires_at_point = dict()
        for child in signal:
            if child.tag != 'wire':
                continue
            start_point = Point2D(float(child.attrib['x1']),
                                  float(child.attrib['y1']))
            end_point = Point2D(float(child.attrib['x2']),
                                float(child.attrib['y2']))
            # if start and end points are the same, i'm not sure what to do
            if start_point == end_point:
                print(child, start_point, end_point)
                print('ERROR: wire start and end points are identical.')
                exit(1)
                continue
            wires.append(child)
            wires_at_point[start_point] = wires_at_point.get(start_point, 0) + 1
            wires_at_point[end_point] = wires_at_point.get(end_point, 0) + 1
        # find corners
        print(wires_at_point)
        corner = []
        for point, connections in wires_at_point.items():
            if connections != 2:
                continue
            print('Corner at %s' % point)


def get_package_points(package_xml):
    """Return the smd and pad points within the given xml package."""
    points = []
    for item in package_xml:
        if item.tag == 'smd':
            x = float(item.attrib['x'])
            y = float(item.attrib['y'])
            points.append(Point2D(x, y))
        elif item.tag == 'pad':
            x = float(item.attrib['x'])
            y = float(item.attrib['y'])
            points.append(Point2D(x, y))
    return points


def read_libraries(filename):
    """Read and return fixed locations of packages within the file."""
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    libraries = dict()
    # store DRU information
    dru = dict()
    for libraries_tag in root.iter('libraries'):
        for item in libraries_tag:
            if item.tag != 'library':
                continue
            library_name = item.attrib['name']
            if library_name not in libraries:
                libraries[library_name] = dict()
            packages = libraries[library_name]
            for package_item in item.iter('package'):
                name = package_item.attrib['name']
                points = get_package_points(package_item)
                packages[name] = points
    return libraries


class Package:
    """A Package holds information on a single component placement."""

    def __init__(self, xml_tag):
        self.origin = Point2D(float(xml_tag.attrib['x']),
                              float(xml_tag.attrib['y']))
        self.library = xml_tag.attrib['library']
        self.footprint = xml_tag.attrib['package']
        if 'rot' in xml_tag.attrib:
            self.mirrored = 'M' in xml_tag.attrib['rot']
            angle = xml_tag.attrib['rot']
            while angle and angle[0].isalpha():
                angle = angle[1:]
            self.rotation = 0.0 if not angle else float(angle)
        else:
            self.mirrored = False
            self.rotation = 0.0

    def __hash__(self):
        return hash((self.library,
                     self.footprint,
                     self.origin,
                     self.rotation,
                     self.mirrored))

    def __repr__(self):
        return ('Package(%s, %s, %s, %s, %s)'
                % (self.library,
                   self.footprint,
                   self.origin,
                   self.rotation,
                   self.mirrored))


def read_placements(filename):
    """
    Return the component placements within the file.

    [('library', 'footprint', x, y, mirror, rotation)]
    """
    # tree = ElementTree.parse(filename)
    # root = tree.getroot()
    placements = [Package(element)
                  for element in
                  ElementTree.parse(filename).getroot().iter('element')]
    return placements


def get_wires_at(layer, point, wires):
    """Return the wires that contain the given point."""
    for wire in wires:
        if wire.attrib['layer'] != layer:
            continue
        start_point = Point2D(float(wire.attrib['x1']),
                              float(wire.attrib['y1']))
        end_point = Point2D(float(wire.attrib['x2']),
                            float(wire.attrib['y2']))
        if start_point == point or end_point == point:
            yield wire


def count_wires_at(layer, point, wires):
    """Return the number of wires that contain the given point."""
    return sum(1 for _ in get_wires_at(layer, point, wires))


def get_wire_endpoints(wire):
    """Return the start and end points of the given wire."""
    start_point = Point2D(float(wire.attrib['x1']),
                          float(wire.attrib['y1']))
    end_point = Point2D(float(wire.attrib['x2']),
                        float(wire.attrib['y2']))
    return start_point, end_point


def get_wire_length(wire):
    """Return the length of the given wire."""
    start_point, end_point = get_wire_endpoints(wire)
    return (start_point - end_point).norm()


def delete_via_teardrops(filename):
    """Delete teardropped vias from the given file."""
    # create script filename
    script_filename = os.path.join(os.path.dirname(filename),
                                   'delete_teardrops.scr')
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    # hold commands to redraw all wires
    commands = []
    commands.append('display none;')
    commands.append('display 1 to 16;')
    commands.append('group all;')
    commands.append('ripup (>0 0);')
    commands.append('display preset_standard;')
    commands.append('set undo_log off;')
    commands.append('set wire_bend 2')
    commands.append('grid mm;')
    # hold wire command to draw by layer
    wires_by_layer = dict()
    # hold all wires that were added during the teardropping step
    teardrop_wires = []
    for signal in root.iter('signal'):
        # read in all wires
        wires = []
        # store vias
        vias = dict()
        # store number of wires at each point and layer
        wires_at_point = dict()
        # store wires
        wires = list()
        for child in signal:
            if child.tag == 'via':
                point = Point2D(float(child.attrib['x']),
                                float(child.attrib['y']))
                assert point not in vias
                vias[point] = child
                continue
            if child.tag != 'wire':
                continue
            start_point = Point2D(float(child.attrib['x1']),
                                  float(child.attrib['y1']))
            end_point = Point2D(float(child.attrib['x2']),
                                float(child.attrib['y2']))
            layer = child.attrib['layer']
            # if start and end points are the same, i'm not sure what to do
            if start_point == end_point:
                print(child, start_point, end_point)
                print('ERROR: wire start and end points are identical.')
                exit(1)
                continue
            wires.append(child)
            key = (start_point, layer)
            wires_at_point[key] = wires_at_point.get(key, 0) + 1
            key = (end_point, layer)
            wires_at_point[key] = wires_at_point.get(key, 0) + 1
        # look for teardropped vias
        for (via_point, layer), wire_count in wires_at_point.items():
            if via_point not in vias:
                continue
            new_teardrop_wires = []
            # a teardrop adds two wires per via per teardrop, so it will
            # have at least 3 wires connected to the via
            if wire_count < 3:
                continue
            print('Teardrop via on signal %s at %s in layer %s with %d '
                  'wires.'
                  % (signal.attrib['name'], via_point, layer, wire_count))
            # look through all wires here
            for via_wire in get_wires_at(layer, via_point, wires):
                # get the end wire point
                p1, p2 = get_wire_endpoints(via_wire)
                if p1 == via_point:
                    other_point = p2
                else:
                    other_point = p1
                # if a teardrop point, should have exactly two wires here
                if count_wires_at(layer, other_point, wires) != 2:
                    continue
                # find the other possible teardrop wire
                both_wires = list(get_wires_at(layer, other_point, wires))
                del both_wires[both_wires.index(via_wire)]
                other_wire = both_wires[0]
                # if this wire isn't curved, this isn't a teardrop wire
                if 'curve' not in other_wire.attrib:
                    continue
                # find the end point
                p1, p2 = get_wire_endpoints(other_wire)
                if p1 == other_point:
                    joining_point = p2
                else:
                    joining_point = p1
                # if this is a joined point should be exactly 4 connections
                if count_wires_at(layer, joining_point, wires) != 4:
                    continue
                # add these wires to the teardrop list
                new_teardrop_wires.append(via_wire)
                new_teardrop_wires.append(other_wire)
            if new_teardrop_wires:
                print('- Found %d teardrop wires to delete'
                      % len(new_teardrop_wires))
                teardrop_wires.extend(new_teardrop_wires)
        for wire in wires:
            layer = wire.attrib['layer']
            if layer not in wires_by_layer:
                wires_by_layer[layer] = []
            if wire in teardrop_wires:
                continue
            # add command to create this wire
            curve = ''
            if 'curve' in wire.attrib:
                curve = wire.attrib['curve']
                if not curve.startswith('-'):
                    curve = '+' + curve
                curve += ' '
            command = ('wire \'%s\' %s (%s %s) %s(%s %s);'
                       % (signal.attrib['name'],
                          wire.attrib['width'],
                          wire.attrib['x1'],
                          wire.attrib['y1'],
                          curve,
                          wire.attrib['x2'],
                          wire.attrib['y2']))
            wires_by_layer[layer].append(command)
    # draw all wires
    for layer in sorted(wires_by_layer.keys()):
        commands.append('layer %s;' % layer)
        commands.extend(wires_by_layer[layer])
    # ratsnest to get rid of airwires
    commands.append('grid last;')
    commands.append('optimize;')
    commands.append('ratsnest;')
    commands.append('set undo_log on;')
    # print('\n' * 3)
    # print('\n'.join(commands))
    with open(script_filename, 'w') as f:
        f.write('\n'.join(commands))
    print('\nScript generated at %s' % (script_filename))


def find_all_via_points(filename):
    """Return a list of all vias in the given board file."""
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    vias = list()
    for signal in root.iter('signal'):
        # read in all wires
        wires = []
        # store number of wires at each point and layer
        wires_at_point = dict()
        # store wires
        wires = list()
        for child in signal:
            if child.tag == 'via':
                point = Point2D(float(child.attrib['x']),
                                float(child.attrib['y']))
                vias.append(point)
    return vias


def get_other_point(wire, start_point):
    """Given the start point and wire, return the end point."""
    p1, p2 = get_wire_endpoints(wire)
    if p1 == start_point:
        return p2
    else:
        assert p2 == start_point
        return p1


def create_signal(p1, a1, p2, a2, signal, width):
    """Return the command to make an arced wire with the given angles."""
    a1.normalize()
    a2.normalize()
    theta = math.acos(a1.dot(a2) / a1.norm() / a2.norm())
    ccw_test = a1.angle() + theta - a2.angle()
    ccw_test = ((ccw_test + math.pi) % (2.0 * math.pi)) - math.pi
    cw_test = a1.angle() - theta - a2.angle()
    cw_test = ((cw_test + math.pi) % (2.0 * math.pi)) - math.pi
    curve = theta
    assert abs(cw_test) < 1e-6 or abs(ccw_test) < 1e-6
    if abs(cw_test) < abs(ccw_test):
        curve = -curve
    command = ('line \'%s\' %s (%s %s) %s (%s %s);'
               % (signal,
                  width,
                  format_mm(p1.x),
                  format_mm(p1.y),
                  format_mm(curve * 180.0 / math.pi, True),
                  format_mm(p2.x),
                  format_mm(p2.y)))
    return command


def format_mm(mm, leading_plus=False):
    """Return a string version of mm."""
    text = '%.8f' % (mm)
    if leading_plus:
        if not text.startswith('-'):
            text = '+' + text
    if '.' in text:
        text = text.rstrip('0').rstrip('.')
    return text


def get_wire_command(wire, signal_name):
    """Return the string command to create the given wire."""
    curve = ''
    if 'curve' in wire.attrib:
        curve = wire.attrib['curve']
        if not curve.startswith('-'):
            curve = '+' + curve
        curve += ' '
    command = ('line \'%s\' %s (%s %s) %s(%s %s);'
               % (signal_name,
                  wire.attrib['width'],
                  wire.attrib['x1'],
                  wire.attrib['y1'],
                  curve,
                  wire.attrib['x2'],
                  wire.attrib['y2']))
    return command


def generate_wire_command(signal_name, width, p1, p2):
    """Return the string command to create the given straight wire."""
    command = ('line \'%s\' %s (%s %s) (%s %s);'
               % (signal_name,
                  width,
                  format_mm(p1.x),
                  format_mm(p1.y),
                  format_mm(p2.x),
                  format_mm(p2.y)))
    return command


def wire_command_sort(command):
    """Return a Point2D of the rightmost point in the given wire command."""
    assert '(' in command
    points = command.split('(')[1:]
    points = [x[:x.index(')')] for x in points]
    points = [Point2D(*[float(x) for x in p.split()]) for p in points]
    return sorted(points)[-1]
    p1, p2 = points
    # print(points)
    # exit(1)
    # p1
    if p1 < p2:
        return p2
    else:
        return p1


def get_angle_deviation(p1, p2, p3):
    """Return the change in angle between p1-p2 and p2-p3."""
    cosine = (p2 - p1).dot(p3 - p2)
    cosine /= (p2 - p1).norm() * (p3 - p2).norm()
    if cosine > 1.0:
        cosine = 1.0
    return math.acos(cosine)


def get_transition_distance(width_mm, angle):
    """Return the target distance to start the transition."""
    # target inner radius of trace as a ratio of trace width
    target_inner_radius_mils = 30
    # maximum distance the trace can move away from the original path (in mils)
    max_trace_deviation_mils = 5
    # target radius based on target inner radius
    radius_1 = target_inner_radius_mils * 0.0254 + width_mm / 2.0
    # target radius based on deviation
    radius_2 = -0.25 * (-2.0 * max_trace_deviation_mils * 0.0254 +
                        width_mm * (math.cos(angle / 2.0) - 1))
    radius_2 /= math.sin(angle / 4.0) ** 2
    length = min(radius_1, radius_2) * math.tan(angle / 2.0)
    length = radius_2 * math.tan(angle / 2.0)
    return length


def create_teardrop(joining_point,
                    tangent,
                    via_point,
                    radius,
                    signal,
                    layer,
                    width):
    """Return the commands for creating a teardrop shape."""
    # tangent should point towards the via
    # (but this is not strictly necessary)
    assert tangent.dot(via_point - joining_point) >= 0.0
    via = via_point - joining_point
    # point = point - via_point
    # get direction normal to tangent
    normal = Point2D(tangent.y, -tangent.x)
    # find d1 and d2
    d1 = 2.0 * (radius - normal.dot(via))
    if d1 == 0.0:
        d1 = 1e100
    else:
        d1 = (radius ** 2 - via.dot(via)) / d1
    d2 = -2.0 * (radius + normal.dot(via))
    if d2 == 0.0:
        d2 = 1e100
    else:
        d2 = (radius ** 2 - via.dot(via)) / d2
    angle = math.atan2(via.y, via.x)
    startAngle1 = angle + math.pi / 2.0
    startAngle2 = angle + math.pi / 2.0
    # startAngle3 = angle + math.pi / 2.0
    if d1 < 0:
        startAngle1 -= math.pi
    if d2 < 0:
        startAngle2 -= math.pi
    endAngle1 = (joining_point - d1 * normal - via_point).angle()
    endAngle2 = (joining_point - d2 * normal - via_point).angle()
    p1 = via_point
    p1 -= radius * Point2D(math.cos(endAngle1), math.sin(endAngle1))
    p2 = via_point
    p2 -= radius * Point2D(math.cos(endAngle2), math.sin(endAngle2))
    # endAngle1 = math.atan2(via_point - d1 * normal)
    # endAngle2 = math.atan2(via_point - d2 * normal)
    commands = []
    # p3 = via_point
    commands.append('line \'%s\' %s (%s %s) %s (%s %s);'
                    % (signal,
                       width,
                       format_mm(joining_point.x),
                       format_mm(joining_point.y),
                       format_mm(endAngle1 - startAngle1, True),
                       format_mm(p1.x),
                       format_mm(p1.y)))
    commands.append('line \'%s\' %s (%s %s) %s (%s %s);'
                    % (signal,
                       width,
                       format_mm(joining_point.x),
                       format_mm(joining_point.y),
                       format_mm(endAngle2 - startAngle2, True),
                       format_mm(p2.x),
                       format_mm(p2.y)))
    return commands
    # endAngle3 = math.atan2(via_point - d3 * normal)
    # if d3 < 0:
    #    startAngle1 -= math.pi
    # get two solutions to the problem
    # p = point + d * normal
    # ||via_point - p|| = r0
    # ndp = normal.dot(point)
    # diff = math.sqrt(-ndp ** 2 + r0 ** 2)
    # d1 = -ndp + diff
    # d2 = -ndp + diff


def snapped_point(point):
    """Return the point snapped to the nearest native resoltion."""
    resolution = 3.125e-6
    return Point2D(math.floor(point.x / resolution + 0.5) * resolution,
                   math.floor(point.y / resolution + 0.5) * resolution)


def round_signals(filename):
    """
    Create rounded signal lines instead of sharp corners.

    First pass:
    * save all straight wires to memory
    * Save curved wires to output command to duplicate
    * Find fixed points due to vias and pads

    Second pass:
    * Find wires/corners to optimize (non-fixed point, 2 wires)
    * For each corner, mark the desired length away from the intersection
    * Find multiple-connected corners (non-fixed point, 3+ wires)

    Third pass:
    * Resolve corners

    """
    # create script filename
    script_filename = os.path.join(os.path.dirname(filename),
                                   'round_signals.scr')
    # if True, will also round corners of 3+ wire junctions
    rounded_junctions = True
    # if True, will enlarge copper traces onto the tPlace plane to help
    # visualize the margin
    # enlarge_to_tplace = True
    # maximum trace width to modify (in mils)
    # min_trace_width_mils = 10
    # target inner radius of trace as a ratio of trace width
    # target_inner_radius_ratio = 1.0
    # if True, polygons will be created for multi-wire junctions
    # else, simple traces will be used
    create_polygons_in_junctions = True
    # if True, create teardrop shapes for vias
    create_teardrop_vias = not True
    # minimum path segment to create
    # min_segment_length_mils = 5
    # target inner radius of trace as a ratio of trace width
    # target_inner_radius_mils = 20
    # maximum distance the trace can move away from the original path (in mils)
    # max_trace_deviation_mils = 3
    # minimum length of path segment to keep
    # min_remaining_segment_mils = 5
    # shortest possible nonzero wire length (native Eagle resolution)
    native_resolution_mm = 3.125e-6
    # if a wire is shorter than this, ignore it
    length_tolerance_mils = 1e-9 / 0.0254
    # parse the XML file
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    # store DRU information
    dru = dict()
    for dru_root in root.iter('designrules'):
        for item in dru_root:
            if item.tag != 'param':
                continue
            dru[item.attrib['name']] = item.attrib['value']
    # convert to mm where possible
    units = dict()
    units['mm'] = 1.0
    units['mic'] = 0.001
    units['mil'] = 0.0254
    units['inch'] = 25.4
    for key, value in dru.items():
        if ' ' in value:
            continue
        for unit in units.keys():
            if value.endswith(unit):
                dru[key] = float(value[:-len(unit)]) * units[unit]
                break
    # hold commands to redraw all wires
    commands = []
    commands.append('set optimizing off;')
    commands.append('display none;')
    commands.append('display 1 to 16;')
    commands.append('group all;')
    commands.append('ripup (>0 0);')
    commands.append('display preset_standard;')
    # commands.append('set undo_log off;')
    commands.append('set wire_bend 2;')
    commands.append('grid mm;')
    # hold straight segments of wire that can be rounded
    # * must be same width
    # * intersections have exactly 2 traces
    # hold all special points which cannot be moved (vias, pads)
    # locked_points = set(find_all_via_points(filename))
    # hold wires by signal
    # wires_by_signal = dict()
    base_locked_points = set()
    # read footprints from board file
    libraries = read_libraries(filename)
    # read part placement in board file
    parts = read_placements(filename)
    # make smd pads locked points
    for part in parts:
        origin = part.origin
        for point in libraries[part.library][part.footprint]:
            point = Point2D(point.x, point.y)
            if part.rotation != 0:
                point.rotate(part.rotation * math.pi / 180.0)
            if part.mirrored:
                point.x = -point.x
            # snap to grid
            base_locked_points.add(origin + point)
    # hold wire commands to draw by layer
    wires_by_layer = dict()
    curved_wire_count = 0
    print('Finding all signal wires.')
    for signal in root.iter('signal'):
        signal_name = signal.attrib['name']
        print('\nIn signal %s:' % (signal_name))
        # if signal_name not in wires_by_signal:
        #    wires_by_signal[signal_name] = []
        # store via points
        via_points = set()
        # store all via tag attributes
        # via[point] = {'x': 1.0, 'y': 1.0, 'drill': x
        vias = dict()
        # store number of wires at each point and layer
        wires_at_point = dict()
        # hold fixed points (vias, end points)
        locked_points = set(base_locked_points)
        # store all straight wires
        wires = list()
        for child in signal:
            if child.tag == 'via':
                p = Point2D(float(child.attrib['x']), float(child.attrib['y']))
                via_points.add(p)
                vias[p] = child.attrib
                locked_points.add(p)
                continue
            if child.tag != 'wire':
                continue
            layer_name = child.attrib['layer']
            start_point = Point2D(float(child.attrib['x1']),
                                  float(child.attrib['y1']))
            end_point = Point2D(float(child.attrib['x2']),
                                float(child.attrib['y2']))
            if 'curve' in child.attrib:
                locked_points.add(start_point)
                locked_points.add(end_point)
                curved_wire_count += 1
                if layer_name not in wires_by_layer:
                    wires_by_layer[layer_name] = []
                command = get_wire_command(child, signal_name)
                wires_by_layer[layer_name].append(command)
                continue
            # wires.append(StraightWire(child, signal.attrib['name']))
            # if start and end points are the same, ignore this wire
            # (not sure if this is possible)
            if start_point == end_point:
                print('WARNING: wire start and end points are identical.')
                locked_points.add(start_point)
                continue
            # if wire widths are different make this a locked point
            wires.append(child)
            # wires_by_signal[signal_name].append(child)
            key = (layer_name, start_point)
            wires_at_point[key] = wires_at_point.get(key, 0) + 1
            key = (layer_name, end_point)
            wires_at_point[key] = wires_at_point.get(key, 0) + 1
        # adjacent_points[(layer, p1)] = [p2, p3, ...]
        adjacent_points = dict()
        for wire in wires:
            p1, p2 = get_wire_endpoints(wire)
            key = (wire.attrib['layer'], p1)
            if key not in adjacent_points:
                adjacent_points[key] = []
            adjacent_points[key].append(p2)
            key = (wire.attrib['layer'], p2)
            if key not in adjacent_points:
                adjacent_points[key] = []
            adjacent_points[key].append(p1)
        # sort adjacent points by angle
        new_adjacent_points = dict()
        for (layer, p), points in adjacent_points.items():
            angles = [(x - p).angle() for x in points]
            items = sorted((x, y) for x, y in zip(angles, points))
            new_points = [x[1] for x in items]
            new_adjacent_points[(layer, p)] = new_points
        adjacent_points = new_adjacent_points
        # create teardrops for vias
        if create_teardrop_vias:
            # loop over all vias
            for (layer, via_point), other_points in adjacent_points.items():
                if via_point not in via_points:
                    continue
                # get diameter of this via
                if 'diameter' in vias[via_point]:
                    diameter = float(vias[via_point]['diameter'])
                else:
                    diameter = 0.0
                # increase based on DRU rules
                drill = float(vias[via_point]['drill'])
                if layer == '1' or layer == '16':
                    dru_diameter = drill * float(dru['rvViaOuter'])
                    dru_diameter = max(dru_diameter, dru['rlMinViaOuter'])
                    dru_diameter = min(dru_diameter, dru['rlMaxViaOuter'])
                else:
                    dru_diameter = drill * float(dru['rvViaInner'])
                    dru_diameter = max(dru_diameter, dru['rlMinViaInner'])
                    dru_diameter = min(dru_diameter, dru['rlMaxViaInner'])
                diameter = max(diameter, dru_diameter)
                wire_width = 0.254
                joining_radius = (diameter - wire_width) / 2.0
                # for each wire coming off of a via, create a teardrop
                for next_point in other_points:
                    length = 0.0
                    last_point = via_point
                    success = False
                    # find the target distance away from the via
                    target_length = 0.030 / 0.0254
                    # length of path so far
                    while True:
                        # get length to next point
                        this_length = (next_point - last_point).norm()
                        if length + this_length > target_length:
                            success = True
                            alpha = (target_length - length) / this_length
                            break
                        # if next point isn't a corner, we can't create the
                        # teardrop TODO
                        temp_point = last_point
                        last_point = next_point
                        length += this_length
                        # find next point
                        next_points = adjacent_points[(layer, last_point)]
                        if len(next_points) != 2:
                            print('Can\'t create teardrop on via')
                            success = False
                            break
                        if next_points[0] == temp_point:
                            next_point = next_points[1]
                        else:
                            assert next_points[1] == temp_point
                            next_point = next_points[0]
                    if not success:
                        continue
                    # create teardrop

                    # def create_teardrop(joining_point,
                    #                    tangent,
                    #                    via_point,
                    #                    radius,
                    #                    signal,
                    #                    layer,
                    #                    width)
                    assert 0 <= alpha <= 1.0
                    teardrop_point = last_point + alpha * (next_point - last_point)
                    tangent = (last_point - next_point).normalize()
                    # normal = Point2D(tangent.y, -tangent.x)
                    print('Creating teardrop from %s to %s'
                          % (teardrop_point, via_point))
                    new_commands = create_teardrop(teardrop_point,
                                                   tangent,
                                                   via_point,
                                                   joining_radius,
                                                   signal_name,
                                                   layer,
                                                   wire_width)
                    if layer not in wires_by_layer:
                        wires_by_layer[layer] = []
                    wires_by_layer[layer].extend(new_commands)
        # hold corner points (2 wires joining)
        corner_points = set()
        # hold wire widths at corner
        # corner_width[(layer, point)] = X
        corner_width = dict()
        # hold junction points (3+ wires intersecting)
        junction_points = set()
        # for each point, hold the distance to the next point we want to fillet
        # rounded_distance[(layer, point)] = X
        rounded_distance = dict()
        # hold segments which will be rounded
        # (layer, point1, point2)
        # rounded_segment = set()
        # hold adjacent points
        # print(locked_points)
        # look through all points with exactly 2 wires and try to simplify
        for (layer, point), wire_count in wires_at_point.items():
            # if this point is fixed, don't modify this intersection
            dist = None
            for p in locked_points:
                this_dist = point.distance_to(p)
                if dist is None or this_dist < dist:
                    dist = this_dist
            #print(point, dist)
            if dist < 5 * native_resolution_mm:
                continue
            if snapped_point(point) in locked_points:
                continue
            # if only one wire, don't modify it
            if wire_count == 1:
                locked_points.add(point)
                continue
            # if more than 2 wires intersect here, mark it as a junction
            if wire_count > 2:
                junction_points.add((layer, point))
                continue
            else:
                corner_points.add((layer, point))
                assert wire_count == 2
            # get the two intersecting wires
            these_wires = list(get_wires_at(layer, point, wires))
            assert len(these_wires) == 2
            [wire_one, wire_two] = these_wires
            # if the two wires are different widths, don't modify them
            if wire_one.attrib['width'] != wire_two.attrib['width']:
                continue
            corner_width[(layer, point)] = wire_one.attrib['width']
            # print('Trying to simplify wires at %s on layer %s'
            #      % (point, layer))
            # print(wires)
            wire_width = float(wire_one.attrib['width'])
            # target_radius = wire_width
            # target_radius *= target_inner_radius_ratio + 0.5
            # target_radius = wire_width + target_inner_radius_mils * 0.0254
            p2 = point
            p1 = get_other_point(wire_one, p2)
            p3 = get_other_point(wire_two, p2)
            assert p3 != p1
            # get the maximum distance based on the path length restriction
            # l1 = get_wire_length(wire_one)
            # l2 = get_wire_length(wire_two)
            # max_distance = l1 if l1 < l2 else l2
            # max_distance -= min_remaining_segment_mils * 0.0254
            # max_distance /= 2.0
            # get unit direction vectors
            a1 = (p2 - p1).normalize()
            a2 = (p3 - p2).normalize()
            cos_theta = a1.dot(a2)
            if cos_theta > 1.0:
                cos_theta = 1.0
            theta = math.acos(cos_theta)
            # skip angles close to 180
            if theta * 180.0 / math.pi > 170.0:
                continue
            # skip very small angles
            if theta * 180.0 / math.pi < 5.0:
                continue
            target_distance = get_transition_distance(wire_width, theta)
            # get maximum radius constraint based on path deviation
            # original_deviation = wire_width / (2.0 * math.cos(theta / 2.0))
            # max_target_radius = max_trace_deviation_mils * 0.0254
            # max_target_radius += original_deviation
            # max_target_radius -= wire_width / 2.0
            # max_target_radius /= (1.0 / math.cos(theta / 2.0) - 1.0)
            # if target_radius > max_target_radius:
            #    target_radius = max_target_radius
            # target_distance = target_radius * math.tan(theta / 2.0)
            # if target_distance > max_distance:
            #    target_distance = max_distance
            key = (layer, p2)
            assert key not in rounded_distance
            rounded_distance[key] = target_distance
            # assert key not in adjacent_point
            # adjacent_point[key] = (p1, p3)
        for (layer, point) in junction_points:
            if not rounded_junctions:
                break
            # get all wires connected to junction
            these_wires = list(get_wires_at(layer, point, wires))
            # if wires are different widths, ignore this point
            wire_widths = set(x.attrib['width'] for x in these_wires)
            if len(wire_widths) > 1:
                continue
            wire_width = float(min(wire_widths))
            assert (layer, point) not in corner_width
            corner_width[(layer, point)] = these_wires[0].attrib['width']
            # sort them by angle from junction to other point
            wire_angles = []
            for wire in these_wires:
                p2 = get_other_point(wire, point)
                wire_angles.append([(p2 - point).angle(), p2])
            these_wires = [(angle, p, wire)
                           for (angle, p), wire in zip(wire_angles,
                                                       these_wires)]
            these_wires.sort()
            last_wire = list(these_wires[0])
            last_wire[0] += 2.0 * math.pi
            these_wires.append(tuple(last_wire))
            # get target radius based on our settings
            target_distance = float('inf')
            # target_radius = wire_width + target_inner_radius_mils * 0.0254
            for i in range(len(these_wires) - 1):
                p1 = these_wires[i][1]
                p3 = these_wires[i + 1][1]
                inner_angle = these_wires[i + 1][0] - these_wires[i][0]
                theta = get_angle_deviation(p1, point, p3)
                if inner_angle * 180.0 / math.pi > 181.0:
                    # print(inner_angle * 180.0 / math.pi, theta)
                    # exit(1)
                    continue
                # skip angles close to 180
                if theta * 180.0 / math.pi > 170.0:
                    continue
                # skip very small angles
                if theta * 180.0 / math.pi < 5.0:
                    continue
                this_distance = get_transition_distance(wire_width, theta)
                target_distance = min(target_distance, this_distance)
                # original_deviation = wire_width / (2.0 * math.cos(theta / 2.0))
                # max_target_radius = max_trace_deviation_mils * 0.0254
                # max_target_radius += original_deviation
                # max_target_radius -= wire_width / 2.0
                # max_target_radius /= (1.0 / math.cos(theta / 2.0) - 1.0)
                # if target_radius > max_target_radius:
                #    target_radius = max_target_radius
            # target_distance = target_radius * math.tan(theta / 2.0)
            key = (layer, point)
            assert key not in rounded_distance
            rounded_distance[key] = target_distance
            # for x in these_wires:
            #    print(x)
            # exit(1)
        # look through each segment and find out length used by transitions
        # committed_length[(layer, p1, p2)] = X
        committed_length = dict()
        for (layer, p1), points in adjacent_points.items():
            if (layer, p1) not in rounded_distance:
                continue
            for p2 in points:
                pa, pb = sorted([p1, p2])
                key = (layer, pa, pb)
                if key not in committed_length:
                    committed_length[key] = 0.0
                committed_length[key] += rounded_distance[(layer, p1)]
        # find scaling factor for each point
        # scaling_factor[(layer, point)] = X
        scaling_factor = dict()
        for (layer, p1, p2), length in committed_length.items():
            segment_length = (p2 - p1).norm()
            if segment_length < length:
                scale = segment_length / length
                key = (layer, p1)
                scaling_factor[key] = min(scaling_factor.get(key, 1.0), scale)
                key = (layer, p2)
                scaling_factor[key] = min(scaling_factor.get(key, 1.0), scale)
        for (layer, point), scale in scaling_factor.items():
            key = (layer, point)
            if key in rounded_distance:
                rounded_distance[key] *= scale
        # if scaling_factor:
        #    print(scaling_factor)
        #    assert False
        # create remaining segment
        for wire in wires:
            p1, p2 = get_wire_endpoints(wire)
            direction = (p2 - p1).normalize()
            layer_name = wire.attrib['layer']
            key = (layer_name, p1)
            pa = p1
            if key in rounded_distance:
                pa += direction * rounded_distance[key]
            pb = p2
            key = (layer_name, p2)
            if key in rounded_distance:
                pb -= direction * rounded_distance[key]
            if layer_name not in wires_by_layer:
                wires_by_layer[layer_name] = []
            wire_length = (pa - pb).norm()
            if wire_length < length_tolerance_mils * 0.0254:
                # print('skipping wire of length %g' % (wire_length))
                continue
            command = generate_wire_command(signal_name,
                                            wire.attrib['width'],
                                            pa,
                                            pb)
            wires_by_layer[layer_name].append(command)
        # create transitions
        for (layer, p2), distance in rounded_distance.items():
            points = adjacent_points[(layer, p2)]
            if len(points) > 2:
                points.append(points[0])
            # create polygon
            if len(points) > 2 and create_polygons_in_junctions:
                # create list of points along with angle for each
                segments = []
                for i in range(len(points) - 1):
                    p1, p3 = points[i], points[i + 1]
                    a1 = (p1 - p2).normalize()
                    a2 = (p3 - p2).normalize()
                    p2a = p2 + a1 * distance
                    p2b = p2 + a2 * distance
                    angle = a2.angle() - a1.angle()
                    if angle < 0.0:
                        angle += 2.0 * math.pi
                    angle = math.pi - angle
                    # angle = math.acos(max(min(a1.dot(a2), 1.0), -1.0))
                    # angle = math.pi - angle
                    segments.append((p2a, -angle * 180.0 / math.pi))
                # create polygon command
                command = ('polygon \'%s\' %s (%s %s)'
                           % (signal_name,
                              corner_width[(layer, p2)],
                              format_mm(segments[-1][0].x),
                              format_mm(segments[-1][0].y)))
                if signal_name == 'N$15':
                    print(points)
                    print(segments)
                    # exit(1)
                for i in range(len(segments)):
                    command += (' %s (%s %s)'
                                % (format_mm(segments[i - 1][1], True),
                                   format_mm(segments[i][0].x),
                                   format_mm(segments[i][0].y)))
                command += ';'
                wires_by_layer[layer].append(command)
                continue
            for i in range(len(points) - 1):
                p1, p3 = points[i], points[i + 1]
                a1 = (p2 - p1).normalize()
                a2 = (p3 - p2).normalize()
                p2a = p2 - a1 * distance
                p2b = p2 + a2 * distance
                new_wire_command = create_signal(p2a, a1, p2b, a2,
                                                 signal_name,
                                                 corner_width[(layer, p2)])
                if layer not in wires_by_layer:
                    wires_by_layer[layer] = []
                wires_by_layer[layer].append(new_wire_command)
        # go through each segment and create transitions and lines
        print('- %d fixed points' % (len(locked_points)))
        print('- %d corners' % (len(corner_points)))
        print('- %d junctions' % (len(junction_points)))
        # print('Needed to scale back %d points' % (len(scaling_factor)))
        # print(rounded_distance)
        # print(adjacent_point)
    # draw all wires
    for layer in sorted(wires_by_layer.keys()):
        commands.append('layer %s;' % layer)
        commands.append('change thermals off;')
        commands.extend(sorted(wires_by_layer[layer], key=wire_command_sort))
        # for x in sorted(wires_by_layer[layer], key=wire_command_sort):
        #    print(wire_command_sort(x))
    # ratsnest to get rid of airwires
    commands.append('grid last;')
    commands.append('optimize;')
    commands.append('set optimizing on;')
    commands.append('ratsnest;')
    # commands.append('set undo_log on;')
    commands.append('group (>0 0);')
    with open(script_filename, 'w') as f:
        f.write('\n'.join(commands))
    print('\nScript generated at %s' % (script_filename))
    # print('\n' * 3)
    # print('\n'.join(commands))


if True:
    libraries = read_libraries(filename)
    parts = read_placements(filename)
    fixed_points = set()
    for part in parts:
        origin = part.origin
        for point in libraries[part.library][part.footprint]:
            point = Point2D(point.x, point.y)
            if part.rotation != 0:
                point.rotate(part.rotation * math.pi / 180.0)
            if part.mirrored:
                point.x = -point.x
            fixed_points.add(origin + point)
    commands = []
    commands.append('display 48;')
    commands.append('change layer 48;')
    commands.append('grid mm;')
    commands.append('change width 0.0254;')
    radius = 0.010 / 0.0254
    for p in fixed_points:
        commands.append(('circle %s %s;' % (p, p + Point2D(radius, 0.0))).replace(',', ''))
    script_filename = os.path.join(os.path.dirname(filename),
                                   'show_fixed_points.scr')
    with open(script_filename, 'w') as f:
        f.write('\n'.join(commands))


round_signals(filename)
# delete_via_teardrops(filename)
