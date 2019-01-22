"""
This script modifies an eagle board file by rounding corners of traces.



Process for creating rounded traces:




"""

import os
import copy
import pdb
import math
from xml.etree import ElementTree

from point2d import Point2D

# filename to modify
L = []
filename = r'C:\Users\tdkostk\Documents\eagle\projects\round_traces\round_traces_test.brd'
filename = r'C:\Users\tdkostk\Documents\eagle\projects\micro_ohmmeter\micro_ohmmeter_rev5.brd'
#filename = r'C:\Users\tdkostk\Documents\eagle\projects\kct-tester\sandia-cable-tester-rev5.brd'
#filename = r'C:\Users\tdkostk\Documents\eagle\projects\kct-tester\sandia-cable-tester-rev5-rounded.brd'
filename = r'C:\Users\tdkostk\Documents\eagle\projects\teardrop_vias\teardrop_test.brd'


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
    target_inner_radius_mils = 20
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
                    width):
    """
    Return the commands for creating a teardrop shape.

    * joining_point is the point on the trace where the teardrop joins.
    * tangent is the tangent direction of the trace at joining_point in the
      direction of the via.
    * via_point is the point of the via
    * radius is the radius of the via where the teardrop should join

    """
    # tangent should point towards the via
    # (but this is not strictly necessary)
    via = via_point - joining_point
    assert tangent.dot(via) >= 0.0
    # point = point - via_point
    # get direction normal to tangent
    normal = Point2D(tangent.y, -tangent.x)
    # d1 = (via.dot(via) - radius ** 2) / (via.dot(normal) - radius)
    # d2 = (via.dot(via) - radius ** 2) / (via.dot(normal) - radius)
    # find d1 and d2
    d1 = 2.0 * (normal.dot(via) - radius)
    if d1 == 0.0:
        d1 = 1e100
    else:
        d1 = (via.dot(via) - radius ** 2) / d1
    d2 = 2.0 * (normal.dot(via) + radius)
    if d2 == 0.0:
        d2 = 1e100
    else:
        d2 = (via.dot(via) - radius ** 2) / d2
    d3 = 2.0 * via.dot(normal)
    if d3 == 0:
        d3 = 1e100
    else:
        d3 = via.dot(via) / d3
    print('\nd1=%s, d2=%s, d3=%s' % (d1, d2, d3))
    print('t=%s, n=%s' % (tangent, normal))
    start_angle_1 = via.angle() + math.pi / 2.0
    start_angle_2 = via.angle() + math.pi / 2.0
    start_angle_3 = via.angle() + math.pi / 2.0
    start_angle_1 = tangent.angle() + math.pi / 2.0
    start_angle_2 = tangent.angle() + math.pi / 2.0
    start_angle_3 = tangent.angle() + math.pi / 2.0
    if d1 < 0:
        start_angle_1 -= math.pi
    if d2 < 0:
        start_angle_2 -= math.pi
    if d3 < 0:
        start_angle_3 -= math.pi
    end_angle_1 = (via - d1 * normal).angle()
    end_angle_2 = (via - d2 * normal).angle()
    end_angle_3 = (via - d3 * normal).angle()
    # find point on radius to join to
    center = joining_point + d1 * normal
    p1 = via_point
    if center.distance_to(p1) < abs(d1):
        p1 = p1 + radius * Point2D(math.cos(end_angle_1), math.sin(end_angle_1))
    else:
        p1 = p1 - radius * Point2D(math.cos(end_angle_1), math.sin(end_angle_1))
    print('d1=%g, dp1=%g' % (d1, center.distance_to(p1)))
    # find second point on radius to join with
    center = joining_point + d2 * normal
    p2 = via_point
    if center.distance_to(p2) < abs(d2):
        p2 = p2 + radius * Point2D(math.cos(end_angle_2), math.sin(end_angle_2))
    else:
        p2 = p2 - radius * Point2D(math.cos(end_angle_2), math.sin(end_angle_2))
    print('d2=%g, dp2=%g' % (d2, center.distance_to(p2)))
    print('ea1=%s, ea2=%s' % (end_angle_1 * 180 / math.pi, end_angle_2 * 180 / math.pi))
    commands = []
    #
    if False:
        p = joining_point + d1 * normal
        commands.append('circle 0.05 (%s %s) (%s %s);'
                        % (format_mm(p.x),
                           format_mm(p.y),
                           format_mm(p.x + 0.3),
                           format_mm(p.y)))
        p = joining_point + d2 * normal
        commands.append('circle 0.05 (%s %s) (%s %s);'
                        % (format_mm(p.x),
                           format_mm(p.y),
                           format_mm(p.x + 0.3),
                           format_mm(p.y)))
    # end_angle_1 = math.atan2(via_point - d1 * normal)
    # end_angle_2 = math.atan2(via_point - d2 * normal)
    # p3 = via_point
    a1 = (end_angle_1 - start_angle_1) * 180 / math.pi
    a1 = ((a1 + 180) % 360) - 180
    a2 = (end_angle_2 - start_angle_2) * 180 / math.pi
    a2 = ((a2 + 180) % 360) - 180
    a3 = (end_angle_3 - start_angle_3) * 180 / math.pi
    a3 = ((a3 + 180) % 360) - 180
    commands.append('line \'%s\' %s (%s %s) %s (%s %s);'
                    % (signal,
                       width,
                       format_mm(joining_point.x),
                       format_mm(joining_point.y),
                       format_mm(a1, True),
                       format_mm(p1.x),
                       format_mm(p1.y)))
    commands.append('line \'%s\' %s (%s %s) %s (%s %s);'
                    % (signal,
                       width,
                       format_mm(joining_point.x),
                       format_mm(joining_point.y),
                       format_mm(a2, True),
                       format_mm(p2.x),
                       format_mm(p2.y)))
    commands.append('line \'%s\' %s (%s %s) %s (%s %s);'
                    % (signal,
                       width,
                       format_mm(joining_point.x),
                       format_mm(joining_point.y),
                       format_mm(a3, True),
                       format_mm(via_point.x),
                       format_mm(via_point.y)))
    return commands
    # endAngle3 = math.atan2(via_point - d3 * normal)
    # if d3 < 0:
    #    start_angle_1 -= math.pi
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
    # if True, will round corners
    round_corners = True
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
    # if True, traces will also be output in junctions to avoid wire stub DRC
    # errors
    create_traces_in_junctions = True
    # if True, create teardrop shapes for vias
    create_teardrop_vias = False
    # if True, will snap points close to the grid points
    snap_to_grid = True
    # tolerance for snapping points, in inces
    snap_tolerance_inch= 1e-6
    # grid spacing in inches
    grid_spacing_inches = 1e-3
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
    commands.append('change thermals off;')
    for layer in sorted(wires_by_layer.keys()):
        commands.append('layer %s;' % layer)
        commands.extend(sorted(wires_by_layer[layer], key=wire_command_sort))
        # for x in sorted(wires_by_layer[layer], key=wire_command_sort):
        #    print(wire_command_sort(x))
    # set view on top layer
    commands.append('change layer 1;')
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


def get_wires_by_signal(filename):
    """Return the wires sorted by signal from the given file."""
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    wire_by_layer = dict()
    # store DRU information
    dru = dict()
    for signals in root.iter('signals'):
        for signal in signals.iter('signal'):
            name = signal.attrib['name']
            for wire in signal.iter('wire'):
                pass


def snap_wires_to_grid(filename, tolerance_inch=1e-6, spacing_inch=1e-3):
    """Snap wires to the grid if they are close."""
    # number of wire endpoints on-grid/snapped/off-grid
    snapped_wires = [0, 0, 0]
    snapped_vias = [0, 0, 0]
    # native resolution in mm
    native_resolution_mm = 3.125e-6
    # max digits per coord
    #decimal_place_resolution = 8
    # hold commands to redraw all wires
    commands = []
    commands.append('set optimizing off;')
    commands.append('display none;')
    commands.append('display 1 to 16 18;')
    commands.append('group all;')
    commands.append('ripup (>0 0);')
    commands.append('display preset_standard;')
    commands.append('set wire_bend 2;')
    commands.append('grid mm;')
    commands.append('change drill %.9f;' % (0.013 * 25.4))
    # search through the XML tree
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    wires_by_layer = dict()
    # store DRU information
    dru = dict()
    for signals in root.iter('signals'):
        for signal in signals.iter('signal'):
            name = signal.attrib['name']
            for wire in signal.iter('wire'):
                on_grid = True
                x1 = float(wire.attrib['x1'])
                y1 = float(wire.attrib['y1'])
                p1 = Point2D(x1, y1)
                x1 /= mm_per_inch
                x1 = math.floor(x1 / spacing_inch + 0.5) * spacing_inch
                x1 *= mm_per_inch
                y1 /= mm_per_inch
                y1 = math.floor(y1 / spacing_inch + 0.5) * spacing_inch
                y1 *= mm_per_inch
                p1_grid = Point2D(x1, y1)
                dist = p1.distance_to(p1_grid)
                if dist / mm_per_inch < tolerance_inch:
                    p1 = p1_grid
                else:
                    on_grid = False
                x2 = float(wire.attrib['x2'])
                y2 = float(wire.attrib['y2'])
                p2 = Point2D(x2, y2)
                x2 /= mm_per_inch
                x2 = math.floor(x2 / spacing_inch + 0.5) * spacing_inch
                x2 *= mm_per_inch
                y2 /= mm_per_inch
                y2 = math.floor(y2 / spacing_inch + 0.5) * spacing_inch
                y2 *= mm_per_inch
                p2_grid = Point2D(x2, y2)
                dist = p2.distance_to(p2_grid)
                if dist / mm_per_inch < tolerance_inch:
                    p2 = p2_grid
                else:
                    on_grid = False
                # snap p1 to the native grid
                x1 = math.floor(p1.x / native_resolution_mm + 0.5) * native_resolution_mm
                x1 = '%.9f' % x1
                y1 = math.floor(p1.y / native_resolution_mm + 0.5) * native_resolution_mm
                y1 = '%.9f' % y1
                # snap p1 to the native grid
                x2 = math.floor(p2.x / native_resolution_mm + 0.5) * native_resolution_mm
                x2 = '%.9f' % x2
                y2 = math.floor(p2.y / native_resolution_mm + 0.5) * native_resolution_mm
                y2 = '%.9f' % y2

                if (float(x1) == float(wire.attrib['x1']) and
                        float(y1) == float(wire.attrib['y1']) and
                        float(x2) == float(wire.attrib['x2']) and
                        float(y2) == float(wire.attrib['y2'])):
                    if on_grid:
                        snapped_wires[0] += 1
                    else:
                        snapped_wires[2] += 1
                else:
                    snapped_wires[1] += 1
                if 'curve' in wire.attrib:
                    curve = wire.attrib['curve'] + ' '
                else:
                    curve = ''
                command = ('wire \'%s\' %s (%s %s) %s(%s %s);'
                           % (signal.attrib['name'],
                              wire.attrib['width'],
                              x1,
                              y1,
                              curve,
                              x2,
                              y2))
                if wire.attrib['layer'] not in wires_by_layer:
                    wires_by_layer[wire.attrib['layer']] = []
                wires_by_layer[wire.attrib['layer']].append(command)
            for via in signal.iter('via'):
                on_grid = True
                x = float(via.attrib['x'])
                y = float(via.attrib['y'])
                p = Point2D(x, y)
                x /= mm_per_inch
                x = math.floor(x / spacing_inch + 0.5) * spacing_inch
                x *= mm_per_inch
                y /= mm_per_inch
                y = math.floor(y / spacing_inch + 0.5) * spacing_inch
                y *= mm_per_inch
                p_grid = Point2D(x, y)
                dist = p.distance_to(p_grid)
                if dist / mm_per_inch < tolerance_inch:
                    p = p_grid
                else:
                    on_grid = False
                x = math.floor(p.x / native_resolution_mm + 0.5) * native_resolution_mm
                x = '%.9f' % x
                y = math.floor(p.x / native_resolution_mm + 0.5) * native_resolution_mm
                y = '%.9f' % y
                if (float(x) == float(via.attrib['x']) and
                        float(y) == float(via.attrib['y'])):
                    if on_grid:
                        snapped_vias[0] += 1
                    else:
                        snapped_vias[2] += 1
                else:
                    snapped_vias[1] += 1
                command = ('via \'%s\' auto round 1-16 (%s %s);'
                           % (signal.attrib['name'],
                              x,
                              y))
                if wire.attrib['layer'] not in wires_by_layer:
                    wires_by_layer[wire.attrib['layer']] = []
                wires_by_layer[wire.attrib['layer']].append(command)
    # print summary
    print('Wire summary:')
    print('- Already on grid: %d' % snapped_wires[0])
    print('- Snapped to grid: %d' % snapped_wires[1])
    print('- Off grid: %d' % snapped_wires[2])
    # print summary
    print('Via summary:')
    print('- Already on grid: %d' % snapped_vias[0])
    print('- Snapped to grid: %d' % snapped_vias[1])
    print('- Off grid: %d' % snapped_vias[2])
    # draw all wires
    for layer in sorted(wires_by_layer.keys()):
        commands.append('layer %s;' % layer)
        commands.extend(sorted(wires_by_layer[layer], key=wire_command_sort))
        # for x in sorted(wires_by_layer[layer], key=wire_command_sort):
        #    print(wire_command_sort(x))
    # set view on top layer
    commands.append('change layer 1;')
    commands.append('grid last;')
    commands.append('optimize;')
    commands.append('set optimizing on;')
    commands.append('ratsnest;')
    commands.append('group (>0 0);')
    # create script filename
    script_filename = os.path.join(os.path.dirname(filename),
                                   'snap_to_grid.scr')
    with open(script_filename, 'w') as f:
        f.write('\n'.join(commands))
    print('\nScript generated at %s' % script_filename)


def read_vias(filename, dru=None):
    """Read in vias from the given file."""
    # search through the XML tree
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    # hold all vias
    vias = []
    # search through the XML tree
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    for signals in root.iter('signals'):
        for signal in signals.iter('signal'):
            name = signal.attrib['name']
            for via in signal.iter('via'):
                vias.append(Via(via, signal_name=signal.attrib['name']))
    # if a DRU is defined, figure out the true of each via
    if dru is not None:
        for via in vias:
            # figure out actual outer diameter
            annular_ring = via.drill * dru.param['rvViaOuter']
            annular_ring = max(annular_ring, dru.param['rlMinViaOuter'])
            annular_ring = min(annular_ring, dru.param['rlMaxViaOuter'])
            dru_diameter = via.drill + 2 * annular_ring
            if hasattr(via, 'diameter'):
                via.outer_diameter = max(via.diameter, dru_diameter)
            else:
                via.outer_diameter = dru_diameter
            # figure out actual inner diameter
            annular_ring = via.drill * dru.param['rvViaInner']
            annular_ring = max(annular_ring, dru.param['rlMinViaInner'])
            annular_ring = min(annular_ring, dru.param['rlMaxViaInner'])
            dru_diameter = via.drill + 2 * annular_ring
            if hasattr(via, 'diameter'):
                via.inner_diameter = max(via.diameter, dru_diameter)
            else:
                via.inner_diameter = dru_diameter
    return vias


def read_design_rules(filename):
    """Return the design result from the given file"""
    # search through the XML tree
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    for dru_xml in root.iter('designrules'):
        return DesignRules(dru_xml)


class Wire:

    def __init__(self, wire, signal_name=None):
        self.p1 = Point2D(float(wire.attrib['x1']), float(wire.attrib['y1']))
        self.p2 = Point2D(float(wire.attrib['x2']), float(wire.attrib['y2']))
        self.signal = signal_name
        self.width = float(wire.attrib['width'])
        self.layer = wire.attrib['layer']
        if 'curve' in wire.attrib:
            self.curve = float(wire.attrib['curve'])
        else:
            self.curve = 0.0

    def get_other_point(self, point):
        """Given one end point, return the other one."""
        if self.p1 == point:
            return self.p2
        assert self.p2 == point
        return self.p1

    def get_curve_points(self):
        """Return the center, radius, and starting angle of the curved Wire."""
        assert self.curve != 0
        distance = self.p1.distance_to(self.p2)
        # sin(theta / 2) == (d / 2) / r
        radius = distance / 2.0 / math.sin(self.curve / 2.0 * math.pi / 180.0)
        radius = abs(radius)
        midpoint = self.p1 + 0.5 * (self.p2 - self.p1)
        normal = (self.p2 - self.p1).rotate(math.pi / 2.0)
        normal.normalize()
        delta = radius ** 2 - midpoint.distance_to(self.p2) ** 2
        if delta < 0.0:
            delta = 0.0
        delta = math.sqrt(delta)
        if self.curve > 0:
            center = midpoint + normal * delta
        else:
            center = midpoint - normal * delta
        start_angle = center.angle_to(self.p1)
        return center, radius, start_angle

    def get_distance_along(self, alpha):
        """
        Return the point alpha percent along the wire.

        return is point, tangent

        """
        assert 0.0 <= alpha <= 1.0
        #if alpha == 0.0:
        #    return self.p1
        #elif alpha == 1.0:
        #    return self.p2
        # process straight line
        if self.curve == 0.0:
            point = self.p1 + alpha * (self.p2 - self.p1)
            tangent = (self.p2 - self.p1).normalize()
            return point, tangent
        #distance = self.p1.distance_to(self.p2)
        # sin(theta / 2) == (d / 2) / r
        #radius = distance / 2.0 / math.sin(self.curve * math.pi / 180.0)
        #midpoint = self.p1 + 0.5 * (self.p2 - self.p1)
        #normal = (self.p2 - self.p1).rotate(math.pi / 2.0)
        #delta = radius ** 2 - midpoint.distance_to(self.p2) ** 2
        #if delta < 0.0:
        #    delta = 0.0
        #delta = math.sqrt(delta)
        #if self.curve > 0:
        #    center = midpoint + normal * delta
        #else:
        #    center = midpoint - normal * delta
        center, radius, start_angle = self.get_curve_points()
        #start_angle = center.angle_to(self.p1)
        angle = start_angle + alpha * self.curve * math.pi / 180
        point = center + radius * Point2D(math.cos(angle), math.sin(angle))
        tangent = Point2D(math.cos(angle + math.pi / 2.0),
                          math.sin(angle + math.pi / 2.0))
        if self.curve < 0:
            tangent = -tangent
        return point, tangent

    def reversed(self):
        """Return an identical wire with the points reversed."""
        wire = copy.copy(self)
        wire.p1, wire.p2 = wire.p2, wire.p1
        wire.curve = -wire.curve
        return wire

    def get_length(self):
        """Return the length of this wire."""
        if self.curve == 0.0:
            return self.p1.distance_to(self.p2)
        _, radius, _ = self.get_curve_points()
        return radius * abs(self.curve) * math.pi / 180.0

    def get_command(self, signal_name=None):
        """Return the command to recreate the wire."""
        command = ('wire \'%s\' %.9f (%.9f %.9f)%s (%.9f %.9f);'
                   % (signal_name,
                      self.width,
                      self.p1.x,
                      self.p1.y,
                      '' if self.curve == 0.0 else '%+.3f ' % self.curve,
                      self.p2.x,
                      self.p2.y))
        return command

    def __repr__(self):
        """Return a string representation."""
        args = []
        for attr in dir(self):
            # skip hidden values/functions
            if attr.startswith('_'):
                continue
            # skip all functions
            if callable(getattr(self, attr)):
                continue
            args.append('%s=%s' % (attr, getattr(self, attr)))
        return '%s(%s)' % (self.__class__.__name__, ', '.join(sorted(args)))


def read_wires(filename):
    """Return all wires in the file."""
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    wires = []
    for signals in root.iter('signals'):
        for signal in signals.iter('signal'):
            for wire in signal.iter('wire'):
                wires.append(Wire(wire, signal_name=signal.attrib['name']))
    return wires


def find_point_on_chain(wire_chain, via_point, target_distance):
    """
    Return the point and tangent to the given wire chain.

    Return value is (point, tangent), (index, alpha).

    where index is the index of the wire this point belongs to and alpha is
    the proportion along that index.

    """
    next_point = via_point
    index = -1
    for wire in wire_chain:
        index += 1
        next_point = wire.get_other_point(next_point)
        if via_point.distance_to(next_point) < target_distance:
            continue
        backwards = (via_point.distance_to(wire.p1) >
                     via_point.distance_to(wire.p2))
        assert not backwards
        low = 0.0
        high = 1.0
        while low != high:
            test = (low + high) / 2.0
            point, tangent = wire.get_distance_along(test)
            if test == low or test == high:
                break
            if via_point.distance_to(point) < target_distance:
                if backwards:
                    high = test
                else:
                    low = test
            else:
                if not backwards:
                    high = test
                else:
                    low = test
        # reverse direction of tangent so it points back to start
        tangent = -tangent
        return (point, tangent), (index, test)
    return None


def create_teardrop_vias(filename):
    """Create a script to convert vias to teardrop vias in the given file."""
    # hold commands to redraw all vias/wires
    commands = []
    commands.append('set optimizing off;')
    commands.append('display none;')
    commands.append('display 1 to 16;')
    commands.append('group all;')
    commands.append('ripup (>0 0);')
    commands.append('display preset_standard;')
    commands.append('set wire_bend 2;')
    commands.append('grid mm;')
    commands.append('change drill %.9f;' % (0.013 * 25.4))
    # hold design rules
    dru = read_design_rules(filename)
    # hold all vias
    vias = read_vias(filename, dru)
    # hold all wires
    wires = read_wires(filename)
    # To figure out wire tolopogy, we first find the number of wires connected
    # to each point/layer.  If there are only two, those wires should be
    # combined to make a wire chain.
    # count number of wires at each layer/point
    wire_count = dict()
    for wire in wires:
        point = (wire.layer, wire.p1)
        wire_count[point] = wire_count.get(point, 0) + 1
        point = (wire.layer, wire.p2)
        wire_count[point] = wire_count.get(point, 0) + 1
    # points with other than 2 wires are fixed points
    end_points = set(key for key, value in wire_count.items() if value != 2)
    # via points are also fixed
    via_points = set(via.origin for via in vias)
    # store points that are midway through a wire chain
    mid_points = set(key
                     for key, value in wire_count.items()
                     if value == 2 and key[1] not in via_points)
    print('Point breakdown:')
    print('- Found %d wires.' % len(wires))
    print('- Found %d wire chain end points.' % len(end_points))
    print('- Found %d via points.' % len(via_points))
    print('- Found %d mid points.' % len(mid_points))
    # store wires by layer/point
    wire_by_point = dict()
    for wire in wires:
        point = (wire.layer, wire.p1)
        wire_by_point[point] = []
        point = (wire.layer, wire.p2)
        wire_by_point[point] = []
    for wire in wires:
        point = (wire.layer, wire.p1)
        wire_by_point[point].append(wire)
        point = (wire.layer, wire.p2)
        wire_by_point[point].append(wire)
    print('wire_by_point:')
    for key, value in wire_by_point.items():
        print('- %s: %s' % (key, value))
    # get map from point to via
    via_at_point = dict()
    for via in vias:
        assert via.origin not in via_at_point
        via_at_point[via.origin] = via
    # figure out wire chains starting at vias
    # wire_chain[(1, Point2D(0, 0))] = [Wire(...), Wire(...)]
    wire_chains = dict()
    for wire in wires:
        for p in [wire.p1, wire.p2]:
            point = (wire.layer, p)
            if point[1] not in via_points:
                continue
            wire_chains[point] = [wire]
            this_chain = wire_chains[point]
            new_point = (wire.layer, wire.get_other_point(point[1]))
            while new_point in mid_points:
                new_wires = wire_by_point[new_point]
                assert len(new_wires) == 2
                if new_wires[0] == this_chain[-1]:
                    this_chain.append(new_wires[1])
                else:
                    assert new_wires[1] == this_chain[-1]
                    this_chain.append(new_wires[0])
                new_point = (new_point[0], this_chain[-1].get_other_point(new_point[1]))
    print('- Found %d wire chains.' % len(wire_chains))
    print('- Wire chain lengths: %s' % [sum(x.get_length() for x in chain)
                                        for chain in wire_chains.values()])
    # change wires within each chain such that points are sorted
    # [wire1, wire2, wire3, ...]
    # via at wire1.p1 with wire1.p2 == wire2.p1, etc.
    for via, chain in wire_chains.items():
        starting_point = via[1]
        new_chain = []
        for wire in chain:
            if wire.p2 == starting_point:
                wire = wire.reversed()
            else:
                assert wire.p1 == starting_point
            new_chain.append(wire)
            # print(via, new_chain, starting_point)
            starting_point = wire.p2
        wire_chains[via] = new_chain
    teardrop_inner_diameter_mm = 0.050 * 25.4
    # for chains which have a via at both ends, delete the second half
    print(via_points)
    print(wire_chains.values())
    for via, chain in wire_chains.items():
        # if there's not a via at the end, keep the entire chain
        if chain[-1].p2 not in via_points:
            continue
        print('Found chain with vias at each end: %s' % chain)
        # get total length of chain
        chain_length = sum(wire.get_length() for wire in chain)
        assert chain_length > 0.0
        # find exact midpoint
        length_so_far = 0.0
        index = 0
        while length_so_far + chain[index].get_length() < chain_length / 2.0:
            length_so_far += chain[index].get_length()
            index += 1
        # find proportion along this segment
        alpha = (chain_length / 2.0 - length_so_far) / chain[index].get_length()
        if alpha < 0.0:
            alpha = 0.0
        elif alpha == 1.0:
            alpha = 1.0
        # delete trail of chain
        if alpha == 0.0:
            chain[:] = chain[:index]
        else:
            chain[:] = chain[:index + 1]
            chain[-1] = copy.copy(chain[-1])
            # TODO: support curved wires
            point, _ = chain[-1].get_distance_along(alpha)
            # chain[-1].p2 = chain[-1].p1 + alpha * (chain[-1].p2 - chain[-1].p1)
            chain[-1].p2 = point
            chain[-1].curve *= alpha
    print('- Wire chain lengths: %s' % [sum(x.get_length() for x in chain)
                                        for chain in wire_chains.values()])
    for x in wire_chains.values():
        print(x)
    # the via diameter must be at least this much more than the wire width
    # in order to create a teardrop via
    tolerance_mm = 0.1
    # hold wire commands to draw by layer
    wires_by_layer = dict()
    for point, chain in wire_chains.items():
        print('Processing chain at point %s: %s' % (point, chain))

        via_point = point[1]
        wire_width = chain[0].width
        via = via_at_point[via_point]
        via_diameter = via.outer_diameter
        total_chain_length = sum(x.get_length() for x in chain)
        print('- Total length: %s' % total_chain_length)
        # if chain is too short, don't do anything
        if total_chain_length < via_diameter / 2.0 + tolerance_mm:
            print('Chain too short: %s' % chain)
            continue
        # can't teardrop vias if the wires are bigger than the via diameter
        if via_diameter <= wire_width + tolerance_mm:
            print('wire at %s too big (%s)' % (via_point, via))
            continue
        r1 = via_diameter / 2.0
        r2 = teardrop_inner_diameter_mm + wire_width / 2.0
        d = math.sqrt((r1 + r2) ** 2 - (r2 + wire_width / 2.0) ** 2)
        result = find_point_on_chain(chain, via_point, d)
        print(result)
        # if chain is not long enough, ignore it
        if result is None:
            print('Chain short but usable: %s' % chain)
            # use end of chain
            print(chain[-1].get_distance_along(1.0))
            result = [list(chain[-1].get_distance_along(1.0)), (len(chain) - 1, 1.0)]
            print(result)
            print(result)
            print(result)
            result[0][1] = -result[0][1]
            #continue
        print(result)
        (junction_point, tangent), (wire_index, alpha) = result
        teardrop_commands = create_teardrop(junction_point,
                                            tangent,
                                            via_point,
                                            (via_diameter - wire_width) / 2.0,
                                            via.signal,
                                            chain[0].width)
        # delete portion of chain between via and junction point
        layer = chain[0].layer
        chain[:] = chain[wire_index:]
        if alpha == 1.0:
            chain[:] = chain[1:]
        elif alpha > 0:
            point, _ = chain[0].get_distance_along(alpha)
            chain[0].curve *= 1.0 - alpha
            chain[0].p1 = point
        for x in teardrop_commands:
            print(x)
        if layer not in wires_by_layer:
            wires_by_layer[layer] = []
        wires_by_layer[layer].extend(teardrop_commands)
        print(via_point, junction_point)
    # add wires in chains
    for point, chain in wire_chains.items():
        if wire.layer not in wires_by_layer:
            wires_by_layer[wire.layer] = []
        for wire in chain:
            wires_by_layer[wire.layer].append(wire.get_command(signal_name=wire.signal))
    for layer in sorted(wires_by_layer.keys()):
        commands.append('layer %s;' % layer)
        commands.extend(sorted(wires_by_layer[layer], key=wire_command_sort))
    # set view on top layer
    commands.append('change layer 1;')
    # ratsnest to get rid of airwires
    commands.append('grid last;')
    commands.append('optimize;')
    commands.append('set optimizing on;')
    commands.append('ratsnest;')
    # commands.append('set undo_log on;')
    commands.append('group (>0 0);')
    # create script filename
    script_filename = os.path.join(os.path.dirname(filename),
                                   'create_teardrops.scr')
    with open(script_filename, 'w') as f:
        f.write('\n'.join(commands))
    print('\nScript generated at %s' % (script_filename))

    #print(wire_chains)
    #if point not in wire_chain:
    #    wire_chain[point] = []

    # figure out wire topology
    # we only care about
    #print(wires)


class DesignRules:

    def __init__(self, dru_xml):
        self.name = dru_xml.attrib['name']
        self.param = dict()
        for param in dru_xml.iter('param'):
            self.param[param.attrib['name']] = param.attrib['value']
        # convert to mm where possible
        units = dict()
        units['mm'] = 1.0
        units['mic'] = 0.001
        units['mil'] = 0.0254
        units['inch'] = 25.4
        for key, value in self.param.items():
            if ' ' in value:
                continue
            for unit in units.keys():
                if value.endswith(unit):
                    self.param[key] = float(value[:-len(unit)]) * units[unit]
                    break
            try:
                self.param[key] = float(value)
            except ValueError:
                pass

    def __repr__(self):
        return 'DesignRules(name=%s, param=%s)' % (self.name, self.param)


class Via:

    def __init__(self, via, signal_name):
        self.signal = signal_name
        self.origin = Point2D(float(via.attrib['x']), float(via.attrib['y']))
        self.drill = float(via.attrib['drill'])
        if 'diameter' in via.attrib:
            self.diameter = float(via.attrib['diameter'])
        self.extent = via.attrib['extent']

    def __repr__(self):
        """Return a string representation."""
        args = []
        for attr in dir(self):
            # skip hidden values/functions
            if attr.startswith('_'):
                continue
            # skip all functions
            if callable(getattr(self, attr)):
                continue
            args.append('%s=%s' % (attr, getattr(self, attr)))
        return '%s(%s)' % (self.__class__.__name__, ', '.join(sorted(args)))




mm_per_inch = 25.4


if False:
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


#snap_wires_to_grid(filename)

round_signals(filename)

create_teardrop_vias(filename)

# delete_via_teardrops(filename)
