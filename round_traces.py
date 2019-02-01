"""
This script modifies an eagle board file by rounding corners of traces.

Usage:
> python round_traces.py [options] [board_filename]

For example:
> python round_traces.py --teardrop /path/to/board.brd
> python round_traces.py --round /path/to/board.brd

"""

# When a file is loaded, first the objects are loaded into memory as follows:
# class Board:
#   board.vias = list of Vias()
#   board.wires = list of Wire()
#   board.dru = DesignRules()
#   board.libraries = dict of library -> LibraryPackage()
#   board.elements = list of Element()

import os
import sys
import copy
import filecmp
import math
from xml.etree import ElementTree
from shutil import copyfile

from point2d import Point2D

####################
# START OF OPTIONS #
####################

# if True, board file will be backed up when creating a script
backup_board_file = True

# if True, will also round corners of 3+ wire junctions
rounded_junctions = True

# if True, polygons will be created for multi-wire junctions
# else, simple traces will be used
polygons_in_junctions = True

# if True, traces will also be output in junctions to avoid wire stub DRC
traces_in_junctions = True

# when creating teardrops, plated diameter must be this much more than
# the signal width or a teardrop will be skipped
teardrop_tolerance_mm = 0.1

# target inner radius of teardrops
teardrop_inner_radius_mm = 0.050 * 25.4
# teardrop_inner_radius_mm = 0.25 * 25.4

# if True, will create polygons for teardrops to avoid unplated regions
create_teardrop_polygons = True

# if True, will also teardrop plated through holes found in packages
create_teardrops_on_pths = True

# when rounding signals, maximum deviation from the original copper
max_trace_deviation_mils = 5

# when rounding signals, targer inner radius
target_inner_radius_mils = 100

# if True, will output more info than normal
verbose = True

# if True, will output way more info than normal
super_verbose = True

##################
# END OF OPTIONS #
##################

# native Eagle positional resolution
native_resolution_mm = 3.125e-6

# mm per inch
mm_per_inch = 25.4

# hold path to Eagle projects
project_directory = r'C:\Users\tdkostk\Documents\eagle\projects'

# project and filename to modify
# board_file = r'\round_traces\round_traces_test.brd'
# board_file = r'\teardrop_vias\teardrop_test.brd'
# board_file = r'\kct-tester\sandia-cable-tester-rev5-rounded.brd'
# board_file = r'\kct-tester\sandia-cable-tester-rev5.brd'
# board_file = r'\sandia_cable_tester\sandia-cable-tester-rev5.brd'
# board_file = r'\sandia_cable_tester\sandia-cable-tester-rev5-round.brd'
# board_file = r'\round_traces\round_traces_test_2.brd'
# board_file = r'\sct-adapters\sct-terminal-block-6.brd'
# board_file = r'\sct-adapters\sct-terminal-block-6-round.brd'
# board_file = r'\sct-adapters\sct-bnc-male.brd'
board_file = r'\micro_ohmmeter\micro_ohmmeter_rev5.brd'

board_file = project_directory + board_file


class SMD:
    """An SMD is a pad for a surface mount component pin."""

    def __init__(self, smd_xml):
        """Initialize using the xml tag."""
        # pad name
        self.name = smd_xml.attrib['name']
        # center of pad
        self.origin = Point2D(float(smd_xml.attrib['x']),
                              float(smd_xml.attrib['y']))
        # width of pad
        self.dx = smd_xml.attrib['dx']
        # height of pad
        self.dy = smd_xml.attrib['dy']
        # layer name
        self.layer = smd_xml.attrib['layer']
        # pad rotation
        if 'rot' in smd_xml.attrib:
            self.rot = float(smd_xml.attrib['rot'][1:])
        else:
            self.rot = 0.0
        # pad roundness
        if 'roundness' in smd_xml.attrib:
            self.roundness = int(smd_xml.attrib['roundness'])
        else:
            self.roundness = 0

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


class PTH:
    """A PTH is a circular plated through hole."""

    def __init__(self,
                 origin=Point2D(0.0, 0.0),
                 drill=0.0,
                 outer_diameter=0.0,
                 inner_diameter=0.0,
                 signal=None):
        # location
        self.origin = origin
        # drill diameter
        self.drill = drill
        # outer diameter on outside layers
        self.outer_diameter = outer_diameter
        # outer diameter on inside layers
        self.inner_diameter = inner_diameter
        # connected signal
        self.signal = signal

    @staticmethod
    def from_via(via):
        """Return a PTH created from the given Via."""
        return PTH(via.origin,
                   via.drill,
                   via.outer_diameter,
                   via.inner_diameter,
                   via.signal)

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


class Pad:
    """A Pad is a plated through hole for a through hole component pin."""

    def __init__(self, pad_xml):
        """Initialize using the xml tag."""
        # pad name
        self.name = pad_xml.attrib['name']
        # center of pad
        self.origin = Point2D(float(pad_xml.attrib['x']),
                              float(pad_xml.attrib['y']))
        # drill diameter
        self.drill = float(pad_xml.attrib['drill'])
        # pad shape
        if 'shape' in pad_xml.attrib:
            self.shape = pad_xml.attrib['shape']
        else:
            self.shape = 'round'
        # pad rotation (for non-circular pads)
        if 'rot' in pad_xml.attrib:
            self.rot = float(pad_xml.attrib['rot'][1:])
        else:
            self.rot = 0.0
        # pad diameter
        if 'diameter' in pad_xml.attrib:
            self.diameter = float(pad_xml.attrib['diameter'])
        else:
            self.diameter = 0.0
        # get first attribute
        if 'first' in pad_xml.attrib:
            self.first = pad_xml.attrib['first']
        else:
            self.first = 'no'

    def get_diameters(self, dru):
        """Return the actual outer and inner diameters given the DRU rules."""
        # figure out actual outer diameter
        annular_ring = self.drill * dru.param['rvViaOuter']
        annular_ring = max(annular_ring, dru.param['rlMinViaOuter'])
        annular_ring = min(annular_ring, dru.param['rlMaxViaOuter'])
        dru_diameter = self.drill + 2 * annular_ring
        outer_diameter = max(self.diameter, dru_diameter)
        # figure out actual inner diameter
        annular_ring = self.drill * dru.param['rvViaInner']
        annular_ring = max(annular_ring, dru.param['rlMinViaInner'])
        annular_ring = min(annular_ring, dru.param['rlMaxViaInner'])
        dru_diameter = self.drill + 2 * annular_ring
        inner_diameter = max(self.diameter, dru_diameter)
        return outer_diameter, inner_diameter

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


class LibraryPackage:
    """A LibraryPackage is a footprint within a library."""

    def __init__(self, package_xml):
        """Initialize using the xml tag."""
        self.pads = []
        self.smds = []
        for item in package_xml:
            if item.tag == 'smd':
                self.smds.append(SMD(item))
            elif item.tag == 'pad':
                self.pads.append(Pad(item))

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


class Polygon:
    """A Polygon holds information about a single polygon within a BRD file."""

    def __init__(self):
        self.signal = None
        self.layer = '1'
        self.width = 0.0
        self.rank = '1'
        self.pour = 'solid'
        self.spacing = 0.0
        self.orphans = 'off'
        self.isolate = 0.0
        self.thermals = 'on'
        self.locked = 'no'
        # self.verticies = [(Point2D(), curve), ...]
        self.vertices = []

    @staticmethod
    def from_xml(polygon_xml):
        """Create a Polygon given an XML tag."""
        assert polygon_xml.tag == 'polygon'
        polygon = Polygon()
        if 'rank' in polygon_xml.attrib:
            polygon.rank = polygon_xml.attrib['rank']
        polygon.width = polygon_xml.attrib['width']
        polygon.layer = polygon_xml.attrib['layer']
        for vertex in polygon_xml.iter('vertex'):
            curve = 0.0
            if 'curve' in vertex.attrib:
                curve = float(vertex.attrib['curve'])
            polygon.vertices.append((Point2D(float(vertex.attrib['x']),
                                             float(vertex.attrib['y'])),
                                     curve))
        print(polygon)
        return polygon

    def get_commands(self):
        """Return the commands needed to generate this polygon."""
        commands = []
        # set rank
        commands.append('change layer %s;' % self.layer)
        commands.append('change pour %s;' % self.pour)
        commands.append('change rank %s;' % self.rank)
        commands.append('change thermals %s;' % self.thermals)
        if self.pour != 'solid':
            commands.append('change spacing %s;'
                            % format_position(self.spacing))
        commands.append('change isolate %s;' % format_position(self.spacing))
        commands.append('change orphans %s;' % self.orphans)
        command = 'polygon'
        if self.signal:
            command += ' \'%s\'' % self.signal
        command += ' %s' % self.width
        for point, curve in self.vertices:
            command += (' (%s %s) %s' %
                        (format_position(point.x),
                         format_position(point.y),
                         format_angle(curve)))
        command += (' (%s %s);' %
                    (format_position(self.vertices[0][0].x),
                     format_position(self.vertices[0][0].y)))
        commands.append(command)
        return commands

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


class Wire:

    def __init__(self,
                 p1=Point2D(),
                 p2=Point2D(),
                 signal=None,
                 width=0.0,
                 layer='1',
                 curve=0.0):
        self.p1 = p1
        self.p2 = p2
        self.signal = signal
        self.width = width
        self.layer = layer
        self.curve = curve

    @staticmethod
    def from_rounded_corner(p1, p2, p3, width, signal, layer):
        """Create and return a rounded wire between the given points."""
        # The distance from p1 to p2 should be the same as from p2 to p3
        # get directions
        a1 = (p2 - p1).normalize()
        a2 = (p3 - p2).normalize()
        cosine = min(1.0, max(-1.0, a1.dot(a2)))
        theta = math.acos(cosine)
        ccw_test = a1.angle() + theta - a2.angle()
        ccw_test = ((ccw_test + math.pi) % (2.0 * math.pi)) - math.pi
        cw_test = a1.angle() - theta - a2.angle()
        cw_test = ((cw_test + math.pi) % (2.0 * math.pi)) - math.pi
        ### print(p1, p2, p3)
        ### print(cosine, theta, a1.angle(), a2.angle(), cw_test, ccw_test)
        assert abs(cw_test) < 1e-6 or abs(ccw_test) < 1e-6
        if abs(cw_test) < abs(ccw_test):
            theta = -theta
        return Wire(p1, p3, signal, width, layer, theta)

    @staticmethod
    def from_xml(wire_xml, signal_name=None):
        """Create and return a wire from the XML tag."""
        assert wire_xml.tag == 'wire'
        wire = Wire()
        wire.p1 = Point2D(float(wire_xml.attrib['x1']),
                          float(wire_xml.attrib['y1']))
        wire.p2 = Point2D(float(wire_xml.attrib['x2']),
                          float(wire_xml.attrib['y2']))
        wire.signal = signal_name
        wire.width = float(wire_xml.attrib['width'])
        wire.layer = wire_xml.attrib['layer']
        if 'curve' in wire_xml.attrib:
            wire.curve = float(wire_xml.attrib['curve']) * math.pi / 180.0
        else:
            wire.curve = 0.0
        return wire

    @staticmethod
    def from_points(p1, p2, width, signal, layer):
        """Create and return a straight wire between the given points."""
        wire = Wire()
        wire.p1 = p1
        wire.p2 = p2
        wire.signal = signal
        wire.width = width
        wire.layer = layer
        wire.curve = 0.0
        return wire

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
        Return the point and tangent alpha percent along the wire.

        return is point, tangent

        """
        assert 0.0 <= alpha <= 1.0
        # process straight line
        if self.curve == 0.0:
            point = self.p1 + alpha * (self.p2 - self.p1)
            tangent = (self.p2 - self.p1).normalize()
            return point, tangent
        #
        center, radius, start_angle = self.get_curve_points()
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

    def get_command(self):
        """Return the command to recreate the wire."""
        if self.curve == 0.0:
            angle = ''
        else:
            angle = ' %s' % format_angle(self.curve * 180.0 / math.pi)
        command = ('wire \'%s\' %s (%s %s)%s (%s %s);'
                   % (self.signal if self.signal else '',
                      format_position(self.width),
                      format_position(self.p1.x),
                      format_position(self.p1.y),
                      angle,
                      format_position(self.p2.x),
                      format_position(self.p2.y)))
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


class Via:

    def __init__(self, via, signal_name):
        self.signal = signal_name
        self.origin = Point2D(float(via.attrib['x']), float(via.attrib['y']))
        self.drill = float(via.attrib['drill'])
        if 'diameter' in via.attrib:
            self.diameter = float(via.attrib['diameter'])
        else:
            self.diameter = 0.0
        self.inner_diameter = self.diameter
        self.outer_diameter = self.diameter
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


class Element:
    """An Element holds information on a single component placement."""

    def __init__(self, xml_tag):
        self.name = xml_tag.attrib['name']
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
        return ('Element(%s, %s, %s, %s, %s)'
                % (self.library,
                   self.footprint,
                   self.origin,
                   self.rotation,
                   self.mirrored))


def simplify_commands(commands):
    """Remove unnecessary commands from the list of commands."""
    parameters = dict()
    # hold unnecessary command lines
    duplicates = []
    for i, command in enumerate(commands):
        command = commands[i].lower()
        if not command.startswith('change '):
            continue
        words = command.split(maxsplit=2)
        if len(words) < 3:
            continue
        if words[1] in parameters and parameters[words[1]] == words[2]:
            duplicates.append(i)
        else:
            parameters[words[1]] = words[2]
    # delete those lines
    for i in reversed(duplicates):
        del commands[i]


class Board:
    """A Board holds information within a single Eagle BRD file."""

    def __init__(self, filename):
        """Initialize with a filename."""
        self.filename = filename
        self.dru = read_design_rules(filename)
        self.libraries = read_libraries(filename)
        self.elements = read_elements(filename)
        self.vias = read_vias(filename, self.dru)
        self.wires = read_wires(filename)
        self.contact_refs = read_contact_refs(filename)
        # Polygons are handled differently than wires since they cannot be
        # easily deleted using command-line options.  Because of this, we store
        # the polygons originally in the board file as well as new polygons to
        # create.
        # polygons in the board file
        self.polygons = read_polygons(filename)
        # new polygons to create
        self.new_polygons = []

    def backup_file(self):
        """Create a board backup file."""
        backup_file(self.filename)

    def generate_script(self, wires=True, polygons=True, vias=False):
        """Generate a script which generates copper in this Board."""
        # hold script commands to run
        commands = []
        commands.append('set optimizing off;')
        commands.append('set wire_bend 2;')
        commands.append('grid mm;')
        ### commands.append('change drill %s;' % format_position(0.013 * 25.4))
        if wires:
            commands.append('display none;')
            commands.append('display 1 to 16;')
            commands.append('group all;')
            commands.append('ripup (>0 0);')
        commands.append('display none;')
        commands.append('display preset_standard;')
        copper_by_layer = dict()
        if wires:
            for wire in self.wires:
                if wire.layer not in copper_by_layer:
                    copper_by_layer[wire.layer] = []
                copper_by_layer[wire.layer].append(wire.get_command())
        # add copper commands
        for layer in sorted(copper_by_layer.keys()):
            commands.append('change layer %s;' % layer)
            commands.extend(copper_by_layer[layer])
        for polygon in self.new_polygons:
            commands.append('change layer %s;' % polygon.layer)
            commands.extend(polygon.get_commands())
        assert not vias
        # create script filename
        script_filename = os.path.join(os.path.dirname(self.filename),
                                       'round_traces.scr')
        if verbose:
            print('- Script generated in board file directory.')
            print('- To run, open board and run "script %s".' % script_filename)
        commands.append('display preset_standard;')
        commands.append('change layer 1;')
        commands.append('grid last;')
        commands.append('optimize;')
        commands.append('set optimizing on;')
        commands.append('ratsnest;')
        # remove unnecessary change commands
        simplify_commands(commands)
        # write script
        with open(script_filename, 'w') as f:
            f.write('\n'.join(commands))
        return script_filename



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
    for libraries_tag in root.iter('libraries'):
        for library_tag in libraries_tag.iter('library'):
            library_name = library_tag.attrib['name']
            assert library_name not in libraries
            libraries[library_name] = dict()
            packages = libraries[library_name]
            for package_tag in library_tag.iter('package'):
                name = package_tag.attrib['name']
                assert name not in packages
                packages[name] = LibraryPackage(package_tag)
    return libraries


def read_elements(filename):
    """
    Return the component placements within the file.

    [('library', 'footprint', x, y, mirror, rotation)]
    """
    # tree = ElementTree.parse(filename)
    # root = tree.getroot()
    placements = [Element(element)
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
    with open(script_filename, 'w') as f:
        f.write('\n'.join(commands))
    print('\nScript generated at %s' % script_filename)


def find_all_via_points(filename):
    """Return a list of all vias in the given board file."""
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    # stora vias
    vias = list()
    for signal in root.iter('signal'):
        for via_xml in signal.iter('via'):
            point = Point2D(float(via_xml.attrib['x']),
                            float(via_xml.attrib['y']))
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
    assert False
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
    command = ('wire \'%s\' %s (%s %s) %s (%s %s);'
               % (signal,
                  width,
                  format_position(p1.x),
                  format_position(p1.y),
                  format_angle(curve * 180.0 / math.pi),
                  format_position(p2.x),
                  format_position(p2.y)))
    return command


def format_position(mm):
    """Return a string position in native resolution."""
    mm = math.floor(mm / native_resolution_mm + 0.5) * native_resolution_mm
    text = '%.9f' % mm
    assert '.' in text
    text = text.rstrip('0').rstrip('.')
    return text


def format_angle(degrees):
    """Return a string version of degrees in native resolution."""
    text = '%+.6f' % degrees
    assert '.' in text
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
    command = ('wire \'%s\' %s (%s %s) %s(%s %s);'
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
    command = ('wire \'%s\' %s (%s %s) (%s %s);'
               % (signal_name,
                  width,
                  format_position(p1.x),
                  format_position(p1.y),
                  format_position(p2.x),
                  format_position(p2.y)))
    return command


def wire_command_sort(command):
    """Return a Point2D of the rightmost point in the given wire command."""
    assert '(' in command
    points = command.split('(')[1:]
    points = [x[:x.index(')')] for x in points]
    points = [Point2D(*[float(x) for x in p.split()]) for p in points]
    return sorted(points)[-1]


def get_angle_deviation(p1, p2, p3):
    """Return the change in angle between p1-p2 and p2-p3."""
    # get unit direction vectors
    a1 = (p2 - p1).normalize()
    a2 = (p3 - p2).normalize()
    cosine = a1.dot(a2)
    # due to roundoff, cosine may be slightly out of bounds
    if cosine > 1.0:
        cosine = 1.0
    elif cosine < -1.0:
        cosine = -1.0
    return math.acos(cosine)


def get_transition_distance(width_mm, p1, p2, p3):
    """Return the target distance to start the transition."""
    # if two of the points are the same, don't round the corner
    if p1 == p3 or p1 == p2 or p2 == p3:
        if verbose:
            print('Info: duplicate points in get_transition_distance()')
        return 0.0
    # get opening angle
    theta = math.pi - get_angle_deviation(p1, p2, p3)
    assert 0 <= theta <= math.pi
    # for angles close to 180, don't round the corner
    #if theta * 180.0 / math.pi > 170.0:
    #    return 0.0
    # skip very small angles
    if theta * 180.0 / math.pi < 5.0:
        print(theta)
        if verbose:
            print('Info: small angle encountered in get_transition_distance()')
        return 0.0
    # target radius based on target inner radius
    radius_1 = target_inner_radius_mils * 0.0254 + width_mm / 2.0
    # target radius based on deviation
    radius_2 = -0.25 * (-2.0 * max_trace_deviation_mils * 0.0254 +
                        width_mm * (math.cos(theta / 2.0) - 1))
    radius_2 /= math.sin(theta / 4.0) ** 2
    length = min(radius_1, radius_2) * math.tan(theta / 2.0)
    assert length > 0.0
    return length


def create_teardrop(joining_point,
                    tangent,
                    via_point,
                    radius,
                    signal,
                    width,
                    polygon_teardrop=False):
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
    # find second point on radius to join with
    center = joining_point + d2 * normal
    p2 = via_point
    if center.distance_to(p2) < abs(d2):
        p2 = p2 + radius * Point2D(math.cos(end_angle_2), math.sin(end_angle_2))
    else:
        p2 = p2 - radius * Point2D(math.cos(end_angle_2), math.sin(end_angle_2))
    commands = []
    # end_angle_1 = math.atan2(via_point - d1 * normal)
    # end_angle_2 = math.atan2(via_point - d2 * normal)
    # p3 = via_point
    a1 = (end_angle_1 - start_angle_1) * 180 / math.pi
    a1 = ((a1 + 180) % 360) - 180
    a2 = (end_angle_2 - start_angle_2) * 180 / math.pi
    a2 = ((a2 + 180) % 360) - 180
    a3 = (end_angle_3 - start_angle_3) * 180 / math.pi
    a3 = ((a3 + 180) % 360) - 180
    # add initial line
    commands.append('wire \'%s\' %s (%s %s) %s (%s %s);'
                    % (signal,
                       width,
                       format_position(joining_point.x),
                       format_position(joining_point.y),
                       format_angle(a3),
                       format_position(via_point.x),
                       format_position(via_point.y)))
    # add polygon or lines
    if polygon_teardrop:
        polygon = ('polygon \'%s\' %s (%s %s) %s (%s %s)'
                   % (signal,
                      width,
                      format_position(joining_point.x),
                      format_position(joining_point.y),
                      format_angle(a1),
                      format_position(p1.x),
                      format_position(p1.y)))
        polygon += (' %s (%s %s) %s (%s %s)'
                    % (format_angle(0),
                       format_position(via_point.x),
                       format_position(via_point.y),
                       format_angle(0),
                       format_position(p2.x),
                       format_position(p2.y)))
        polygon += (' %s (%s %s);'
                    % (format_angle(-a2),
                       format_position(joining_point.x),
                       format_position(joining_point.y)))
        commands.append(polygon)
    else:
        commands.append('wire \'%s\' %s (%s %s) %s (%s %s);'
                        % (signal,
                           width,
                           format_position(joining_point.x),
                           format_position(joining_point.y),
                           format_angle(a1),
                           format_position(p1.x),
                           format_position(p1.y)))
        commands.append('wire \'%s\' %s (%s %s) %s (%s %s);'
                        % (signal,
                           width,
                           format_position(joining_point.x),
                           format_position(joining_point.y),
                           format_angle(a2),
                           format_position(p2.x),
                           format_position(p2.y)))
    return commands


def snapped_point(point):
    """Return the point snapped to the nearest native resolution."""
    return Point2D(float(format_position(point.x)),
                   float(format_position(point.y)))


def backup_file(filename):
    """Create a backup of a board file."""
    # ensure file exists
    assert os.path.isfile(filename)
    # get base board directory
    dir_name = os.path.dirname(filename)
    # get base board filename
    base_name = os.path.basename(filename)
    if not base_name.lower().endswith('.brd'):
        print('WARNING: Board file doesn\'t have expected BRD extension.')
    # max backup files
    maximum_backup_file_count = 99
    # find the next available backup name
    for i in range(1, maximum_backup_file_count + 1):
        backup_name = 'backup_%s.%d' % (base_name, i)
        backup_path = os.path.join(dir_name, backup_name)
        # if file is the same, no need to backup again
        if not os.path.isfile(backup_path):
            copyfile(filename, backup_path)
            if not os.path.isfile(backup_path):
                print('ERROR: Could not create board file backup.')
            print('- Backup board file "%s" created.' % backup_name)
            return
        if filecmp.cmp(filename, backup_path):
            print('- Backup board file already exists.')
            return
    print('WARNING: Backup file limit exceeded.  No backup created.')


def obsolete_round_signals(filename):
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
    base_board_filename = os.path.basename(filename)
    print('\n\nRounding traces in board file "%s".' % base_board_filename)
    # create script filename
    script_filename = os.path.join(os.path.dirname(filename),
                                   'round_signals.scr')
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
    # hold points which should not be moved
    base_locked_points = set()
    # read footprints from board file
    libraries = read_libraries(filename)
    # read part placement in board file
    parts = read_elements(filename)
    # make smd pads locked points
    for part in parts:
        origin = part.origin
        footprint = libraries[part.library][part.footprint]
        for point in [*footprint.pads, *footprint.smds]:
            point = Point2D(point.origin.x, point.origin.y)
            if part.rotation != 0:
                point.rotate(part.rotation * math.pi / 180.0)
            if part.mirrored:
                point.x = -point.x
            # snap to grid
            base_locked_points.add(origin + point)
    # hold number of corners rounded
    rounded_corner_count = 0
    # hold number of junctions rounded
    rounded_junction_count = 0
    # hold wire commands to draw by layer
    wires_by_layer = dict()
    curved_wire_count = 0
    for signal in root.iter('signal'):
        signal_name = signal.attrib['name']
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
        # look through all points with exactly 2 wires and try to simplify
        for (layer, point), wire_count in wires_at_point.items():
            # if this point is fixed, don't modify this intersection
            dist = None
            for p in locked_points:
                this_dist = point.distance_to(p)
                if dist is None or this_dist < dist:
                    dist = this_dist
            if dist is not None and dist < 5 * native_resolution_mm:
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
            wire_width = float(wire_one.attrib['width'])
            p2 = point
            p1 = get_other_point(wire_one, p2)
            p3 = get_other_point(wire_two, p2)
            assert p3 != p1
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
                    continue
                # skip angles close to 180
                if theta * 180.0 / math.pi > 170.0:
                    continue
                # skip very small angles
                if theta * 180.0 / math.pi < 5.0:
                    continue
                this_distance = get_transition_distance(wire_width, theta)
                target_distance = min(target_distance, this_distance)
            # target_distance = target_radius * math.tan(theta / 2.0)
            key = (layer, point)
            assert key not in rounded_distance
            rounded_distance[key] = target_distance
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
            if len(points) > 2 and polygons_in_junctions:
                rounded_junction_count += 1
                # create list of points along with angle for each
                segments = []
                for i in range(len(points) - 1):
                    p1, p3 = points[i], points[i + 1]
                    a1 = (p1 - p2).normalize()
                    a2 = (p3 - p2).normalize()
                    p2a = p2 + a1 * distance
                    angle = a2.angle() - a1.angle()
                    if angle < 0.0:
                        angle += 2.0 * math.pi
                    angle = math.pi - angle
                    segments.append((p2a, -angle * 180.0 / math.pi))
                # create polygon command
                if polygons_in_junctions:
                    command = ('polygon \'%s\' %s (%s %s)'
                               % (signal_name,
                                  corner_width[(layer, p2)],
                                  format_position(segments[-1][0].x),
                                  format_position(segments[-1][0].y)))
                    for i in range(len(segments)):
                        command += (' %s (%s %s)'
                                    % (format_angle(segments[i - 1][1]),
                                       format_position(segments[i][0].x),
                                       format_position(segments[i][0].y)))
                    command += ';'
                    wires_by_layer[layer].append(command)
                if (not polygons_in_junctions or
                        traces_in_junctions):
                    for i in range(len(segments) - 1):
                        command = ('wire \'%s\' %s (%s %s) %s (%s %s);'
                                   % (signal_name,
                                      corner_width[(layer, p2)],
                                      format_position(segments[i][0].x),
                                      format_position(segments[i][0].y),
                                      format_angle(segments[i][1]),
                                      format_position(segments[i + 1][0].x),
                                      format_position(segments[i + 1][0].y)))
                        wires_by_layer[layer].append(command)
                continue
            rounded_corner_count += 1
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
    print('- Rounded %d corners and %d junctions.'
          % (rounded_corner_count, rounded_junction_count))
    # draw all wires
    commands.append('change thermals off;')
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
    # backup board
    if backup_board_file:
        backup_file(filename)
    # create script
    with open(script_filename, 'w') as f:
        f.write('\n'.join(commands))
    print('- Script generated in board file directory.')
    print('- To run, open board and run "script %s".' % script_filename)


def round_signals(board):
    """
    Round signals within the given Board object.

    """
    if verbose:
        print('\nRounding signals in board')
    # hold positions of PTHs which should not be moved for each signal
    # locked_pths[signal_name] = {Point2D(), ...}
    locked_pths = {signal: {x.origin for x in pths}
                   for signal, pths in read_pths_by_signal(board).items()}
    if verbose:
        print('- Found %d PTHs for %d signals'
              % (sum(len(x) for x in locked_pths.values()),
                 len(locked_pths)))
    # print(locked_pths)
    # hold points which should not be moved for each signal
    # locked_points[signal_name] = {(layer_name, Point2D), ...}
    locked_points = read_smds_by_signal(board)
    if verbose:
        print('- Found %d SMD pads for %d signals'
              % (sum(len(x) for x in locked_points.values()),
                 len(locked_points)))
    # print(locked_points)
    # hold number of corners rounded
    rounded_corner_count = 0
    # hold number of junctions rounded
    rounded_junction_count = 0
    # number of wires which were already curved
    curved_wire_count = 0
    # hold list of all signals to pass through
    signal_names = set(wire.signal for wire in board.wires)
    # sort wires by signal
    wires_by_signal = dict()
    for signal in signal_names:
        wires_by_signal[signal] = []
    for wire in board.wires:
        wires_by_signal[wire.signal].append(wire)
    # hold list of wires to draw to replace all wires corrently on board
    new_wires = []
    # hold list of new polygons
    new_polygons = []
    # loop through each signal
    for signal_name in sorted(signal_names):
        if super_verbose:
            print('- Rounding wires on signal "%s"' % signal_name)
        # alias locked pths and points
        these_locked_points = locked_points.get(signal_name, {})
        these_locked_pths = locked_pths.get(signal_name, {})
        # alias wires for this signal
        these_wires = wires_by_signal[signal_name]
        if super_verbose:
            print('  - Found %d wires' % len(these_wires))
        # store wires at each (layer, point)
        wires_at_point = dict()
        for wire in these_wires:
            for point in [wire.p1, wire.p2]:
                this_point = (wire.layer, point)
                if this_point not in wires_at_point:
                    wires_at_point[this_point] = []
                wires_at_point[this_point].append(wire)
        # signal_name = signal.attrib['name']
        # store via points
        ### via_points = set()
        # store all via tag attributes
        # via[point] = {'x': 1.0, 'y': 1.0, 'drill': x
        ### vias = dict()
        # hold fixed points (vias, end points)
        ### locked_points = set(base_locked_points)
        # for each wire point, hold the adjacent points
        # adjacent_points[(layer, p1)] = [p2, p3, ...]
        adjacent_points = dict()
        for wire in these_wires:
            if wire.p1 == wire.p2:
                continue
            key = (wire.layer, wire.p1)
            if key not in adjacent_points:
                adjacent_points[key] = []
            adjacent_points[key].append(wire.p2)
            key = (wire.layer, wire.p2)
            if key not in adjacent_points:
                adjacent_points[key] = []
            adjacent_points[key].append(wire.p1)
        # sort adjacent points by angle
        new_adjacent_points = dict()
        for (layer, p), points in adjacent_points.items():
            angles = [(x - p).angle() for x in points]
            items = sorted((x, y) for x, y in zip(angles, points))
            new_points = [x[1] for x in items]
            new_adjacent_points[(layer, p)] = new_points
        adjacent_points = new_adjacent_points
        # hold wire widths at corner
        # corner_width[(layer, point)] = X
        corner_width = dict()
        # for wires of different widths, make the corner a locked point
        # for wires which have the same start and end point, make it locked
        # for curved wires, make both ends a locked point
        for wire in these_wires:
            if wire.p1 == wire.p2:
                locked_points[signal_name].add((wire.layer, wire.p1))
                continue
            if wire.curve != 0.0:
                locked_points[signal_name].add((wire.layer, wire.p1))
                locked_points[signal_name].add((wire.layer, wire.p2))
                continue
            for point in [wire.p1, wire.p2]:
                key = (wire.layer, point)
                if key not in corner_width:
                    corner_width[key] = wire.width
                else:
                    if corner_width[key] != wire.width:
                        locked_points[signal_name].add(key)
        # hold corner points (2 wires joining)
        corner_points = set()
        # hold junction points (3+ wires intersecting)
        junction_points = set()
        # look through all points to find corner and junction points
        if super_verbose:
            print('  - Found %d locked PTHs: %s'
                  % (len(these_locked_pths), these_locked_pths))
            print('  - Found %d locked points: %s'
                  % (len(these_locked_points), these_locked_points))
        for key in wires_at_point.keys():
            layer, point = key
            wire_count = len(wires_at_point[key])
            if super_verbose:
                print('  - Looking at point %s on layer %s with %d neighbors'
                      % (point, layer, wire_count))
            # if only one wire, don't modify it
            if wire_count == 1:
                # I don't think I need to lock it
                # locked_points[signal_name].add(point)
                if super_verbose:
                    print('    - Only one neighbor, skipping')
                continue
            # check for exact match to locked point/pth
            if point in these_locked_pths:
                if super_verbose:
                    print('    - Point is a locked PTH')
                continue
            if (layer, point) in these_locked_points:
                if super_verbose:
                    print('    - Point is a locked point')
                continue
            # get distance to nearest locked point or pth
            dist = float('inf')
            for locked_layer, locked_point in these_locked_points:
                if locked_layer != layer:
                    continue
                this_dist = point.distance_to(locked_point)
                if this_dist < dist:
                    dist = this_dist
            for locked_point in these_locked_pths:
                this_dist = point.distance_to(locked_point)
                if this_dist < dist:
                    dist = this_dist
            # if distance to nearest locked point/pth is small, skip this point
            if dist < 5 * native_resolution_mm:
                if super_verbose:
                    if dist == 0.0:
                        print('    - Point is a locked PTH/point')
                    else:
                        print('    - Point is close to a locked PTH/point')
                continue
            # if more than 2 wires intersect here, mark it as a junction
            if wire_count > 2:
                if super_verbose:
                    print('    - Point marked as a junction')
                junction_points.add(key)
                continue
            # mark it as a corner point to be rounded later
            if super_verbose:
                print('    - Point marked as a corner')
            assert wire_count == 2
            corner_points.add(key)
        # for each point, hold the distance to the next point we want to fillet
        # rounded_distance[(layer, point)] = X
        rounded_distance = dict()
        # get distance to round each corner point
        for (layer, point) in corner_points:
            # get the two intersecting wires
            local_wires = wires_at_point[(layer, point)]
            assert len(local_wires) == 2
            [wire_one, wire_two] = local_wires
            # wires should have same width
            assert wire_one.width == wire_two.width
            assert wire_one.width == corner_width[(layer, point)]
            p2 = point
            p1 = wire_one.get_other_point(p2)
            p3 = wire_two.get_other_point(p2)
            assert p1 != p3
            distance = get_transition_distance(wire_one.width, p1, p2, p3)
            key = (layer, p2)
            assert key not in rounded_distance
            rounded_distance[key] = distance
        # get distance to round each junction point
        for (layer, point) in junction_points:
            if not rounded_junctions:
                break
            # get all wires connected to junction
            local_wires = wires_at_point[(layer, point)]
            # if wires are different widths, ignore this point
            wire_widths = set(x.width for x in local_wires)
            if len(wire_widths) > 1:
                continue
            wire_width = float(min(wire_widths))
            assert (layer, point) in corner_width
            assert corner_width[(layer, point)] == local_wires[0].width
            ### corner_width[(layer, point)] = local_wires[0].attrib['width']
            # sort wires by angle from junction to the other point
            wire_angles = []
            for wire in local_wires:
                p2 = wire.get_other_point(point)
                wire_angles.append([(p2 - point).angle(), p2])
            wire_spokes = [(angle, p, wire)
                           for (angle, p), wire in zip(wire_angles,
                                                       local_wires)]
            # TODO: possible issue with sorting a Wire not implemented
            # It shouldn't be needed, though.  I ran into this issue when
            # running this on a file which had been teardropped.
            wire_spokes.sort()
            last_wire = list(wire_spokes[0])
            last_wire[0] += 2.0 * math.pi
            wire_spokes.append(tuple(last_wire))
            # get target radius based on our settings
            target_distance = float('inf')
            # target_radius = wire_width + target_inner_radius_mils * 0.0254
            for i in range(len(wire_spokes) - 1):
                p1 = wire_spokes[i][1]
                p3 = wire_spokes[i + 1][1]
                distance = get_transition_distance(wire_width, p1, point, p3)
                ### print(distance, p1, point, p3)
                # assert distance < 1000
                assert distance > 0.0
                target_distance = min(target_distance, distance)
            # target_distance = target_radius * math.tan(theta / 2.0)
            key = (layer, point)
            assert key not in rounded_distance
            rounded_distance[key] = target_distance
        if False and signal_name == 'N$1':
            print('\n'.join('%s' % x for x in these_wires))
            print('these_locked_points=%s' % these_locked_points)
            print('these_locked_pths=%s' % these_locked_pths)
            print('adjacent_points=%s' % adjacent_points)
            print('corner_points=%s' % corner_points)
            print('junction_points=%s' % junction_points)
            print('rounded_distance=%s' % rounded_distance)
            exit(1)
        #if junction_points:
        #    exit(1)
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
        # scale down rounding distance where necessary
        for key, scale in scaling_factor.items():
            if scale == 1.0:
                continue
            if key in rounded_distance:
                rounded_distance[key] *= scale
        # for very small radii, don't round the wire at all
        for key in sorted(rounded_distance.keys()):
            if rounded_distance[key] < native_resolution_mm * 10:
                del rounded_distance[key]
        # create remaining part of existing straight segments
        for wire in these_wires:
            # if not modified at all, just copy the wire
            if ((wire.layer, wire.p1) not in rounded_distance and
                    (wire.layer, wire.p2) not in rounded_distance):
                new_wires.append(wire)
                continue
            direction = (wire.p2 - wire.p1).normalize()
            ### layer_name = wire.attrib['layer']
            key = (wire.layer, wire.p1)
            pa = copy.copy(wire.p1)
            if key in rounded_distance:
                pa += direction * rounded_distance[key]
            pb = copy.copy(wire.p2)
            key = (wire.layer, wire.p2)
            if key in rounded_distance:
                pb -= direction * rounded_distance[key]
            ### if layer_name not in wires_by_layer:
            ###     wires_by_layer[layer_name] = []
            # if entire wire is rounded, the middle doesn't need drawn
            wire_length = (pa - pb).norm()
            if wire_length < native_resolution_mm:
                continue
            # generate a new wire
            ### command = generate_wire_command(signal_name,
            ###                                  wire.attrib['width'],
            ###                                pa,
            ###                                pb)
            new_wires.append(Wire.from_points(pa,
                                              pb,
                                              wire.width,
                                              wire.signal,
                                              wire.layer))
        # create transitions
        ### print(sorted(rounded_distance.keys()))
        ### print(sorted(adjacent_points.keys()))
        # be present in adjacent_points unless the latter is deep copied
        # adjacent_points = copy.deepcopy(adjacent_points)
        # corner_width = copy.deepcopy(corner_width)
        # rounded_distance = copy.deepcopy(rounded_distance)
        # adjacent_points = {key: value for key, value in adjacent_points.items()}
        for (layer, p2), distance in rounded_distance.items():
            # print(sorted(adjacent_points.keys()).index((layer, p2)))
            # assert (layer, p2) in adjacent_points
            points = adjacent_points[(layer, p2)]
            # it's a junction point, so duplicate the starting point at the end
            if len(points) > 2:
                points.append(points[0])
            # create polygon
            if len(points) > 2 and polygons_in_junctions:
                rounded_junction_count += 1
                # create list of points along with angle for each
                segments = []
                for i in range(len(points) - 1):
                    p1, p3 = points[i], points[i + 1]
                    a1 = (p1 - p2).normalize()
                    a2 = (p3 - p2).normalize()
                    p2a = p2 + a1 * distance
                    angle = a2.angle() - a1.angle()
                    if angle < 0.0:
                        angle += 2.0 * math.pi
                    angle = math.pi - angle
                    segments.append((p2a, -angle * 180.0 / math.pi))
                # create polygon command
                if polygons_in_junctions:
                    # command = ('polygon \'%s\' %s (%s %s)'
                    #            % (signal_name,
                    #               corner_width[(layer, p2)],
                    #               format_position(segments[-1][0].x),
                    #               format_position(segments[-1][0].y)))
                    # for i in range(len(segments)):
                    #     command += (' %s (%s %s)'
                    #                 % (format_angle(segments[i - 1][1]),
                    #                    format_position(segments[i][0].x),
                    #                    format_position(segments[i][0].y)))
                    # command += ';'
                    ### wires_by_layer[layer].append(command)
                    polygon = Polygon()
                    polygon.thermals = 'off'
                    polygon.width = corner_width[(layer, p2)]
                    polygon.signal = signal_name
                    polygon.layer = layer
                    polygon.vertices.extend(segments)
                    new_polygons.append(polygon)
                if traces_in_junctions or not polygons_in_junctions:
                    for i in range(len(segments) - 1):
                        break
                        new_wires.append(Wire(p1=segments[i][0],
                                              p2=segments[i + 1][0],
                                              signal=signal_name,
                                              width=corner_width[(layer, p2)],
                                              layer=layer,
                                              curve=segments[i][1]))
                        # command = ('wire \'%s\' %s (%s %s) %s (%s %s);'
                        #            % (signal_name,
                        #               corner_width[(layer, p2)],
                        #               format_position(segments[i][0].x),
                        #               format_position(segments[i][0].y),
                        #               format_angle(segments[i][1]),
                        #               format_position(segments[i + 1][0].x),
                        #               format_position(segments[i + 1][0].y)))
                        ### wires_by_layer[layer].append(command)
                continue
            rounded_corner_count += 1
            for i in range(len(points) - 1):
                p1, p3 = points[i], points[i + 1]
                a1 = (p2 - p1).normalize()
                a2 = (p3 - p2).normalize()
                p2a = p2 - a1 * distance
                p2b = p2 + a2 * distance
                new_wire = Wire.from_rounded_corner(
                    p2a,
                    p2,
                    p2b,
                    width=corner_width[(layer, p2)],
                    signal=signal_name,
                    layer=layer)
                new_wires.append(new_wire)
                #new_wire_command = create_signal(p2a, a1, p2b, a2,
                #                                 signal_name,
                #                                 corner_width[(layer, p2)])
                ###if layer not in wires_by_layer:
                ###    wires_by_layer[layer] = []
                ### wires_by_layer[layer].append(new_wire_command)
        # go through each segment and create transitions and lines
    # replace wires with new wires
    if verbose:
        print('- Rounded %d corners and %d junctions.'
              % (rounded_corner_count, rounded_junction_count))
        print('- Replaced %d original wires with %d new wires'
              % (len(board.wires), len(new_wires)))
        print('- Augmented %d original polygons with %d new polygons'
              % (len(board.polygons), len(new_polygons)))
    # replace wires with new wires
    board.wires = new_wires
    ### TODO: remove
    print('\n'.join('%s' % x for x in board.wires if x.signal == 'N$1'))
    # add created polygons
    board.new_polygons.extend(new_polygons)
    # draw all wires
    # commands.append('change thermals off;')
    # for layer in sorted(wires_by_layer.keys()):
    #     commands.append('layer %s;' % layer)
    #     commands.extend(sorted(wires_by_layer[layer], key=wire_command_sort))
    # # set view on top layer
    # commands.append('change layer 1;')
    # # ratsnest to get rid of airwires
    # commands.append('grid last;')
    # commands.append('optimize;')
    # commands.append('set optimizing on;')
    # commands.append('ratsnest;')
    # # commands.append('set undo_log on;')
    # commands.append('group (>0 0);')
    # # backup board
    # if backup_board_file:
    #     backup_file(filename)
    # # create script
    # with open(script_filename, 'w') as f:
    #     f.write('\n'.join(commands))
    # print('- Script generated in board file directory.')
    # print('- To run, open board and run "script %s".' % script_filename)


def snap_wires_to_grid(filename, tolerance_inch=1e-6, spacing_inch=1e-3):
    """Snap wires to the grid if they are close."""
    print('WARNING: this function was never finished.')
    # number of wire endpoints on-grid/snapped/off-grid
    snapped_wires = [0, 0, 0]
    snapped_vias = [0, 0, 0]
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
    commands.append('change drill %s;' % format_position(0.013 * 25.4))
    # search through the XML tree
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    wires_by_layer = dict()
    for signals in root.iter('signals'):
        for signal in signals.iter('signal'):
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
                x1 = format_position(p1.x)
                y1 = format_position(p1.y)
                # snap p1 to the native grid
                x2 = format_position(p2.x)
                y2 = format_position(p2.y)
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
                x = format_position(p.x)
                y = format_position(p.x)
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
                if '1' not in wires_by_layer:
                    wires_by_layer['1'] = []
                wires_by_layer['1'].append(command)
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
    """Read and return vias from the given file."""
    # hold all vias
    vias = []
    # search through the XML tree
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    for signals in root.iter('signals'):
        for signal in signals.iter('signal'):
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


def read_contact_refs(filename):
    """Read and return contact references from the given file."""
    # if LED1 pad A is connected to VCC
    # --> contact_ref[('LED1', 'A')] = 'VCC'
    contact_ref = dict()
    # search through the XML tree
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    for signal in root.iter('signal'):
        signal_name = signal.attrib['name']
        for ref in signal.iter('contactref'):
            key = (ref.attrib['element'], ref.attrib['pad'])
            assert key not in contact_ref
            contact_ref[key] = signal_name
    return contact_ref


def read_pths_by_signal(board):
    """Read and return pths by signal for the given board."""
    pths = dict()
    # add vias for each signal
    for via in board.vias:
        if via.signal not in pths:
            pths[via.signal] = []
        pths[via.signal].append(PTH.from_via(via))
    # look at each package and add PTHs found in library
    for element in board.elements:
        # store element origin
        origin = element.origin
        # alias this footprint
        footprint = board.libraries[element.library][element.footprint]
        # add an entry for each pad in the footprint
        for pad in footprint.pads:
            point = copy.copy(pad.origin)
            if element.rotation != 0:
                point.rotate(element.rotation * math.pi / 180.0)
            if element.mirrored:
                point.x = -point.x
            outer, inner = pad.get_diameters(board.dru)
            if (element.name, pad.name) in board.contact_refs:
                signal = board.contact_refs[(element.name, pad.name)]
                if signal not in pths:
                    pths[signal] = []
                this_pth = PTH(origin + point, pad.drill, outer, inner, signal)
                pths[signal].append(this_pth)
    return pths


def read_smds_by_signal(board):
    """Read and return SMD pads by signal for the given board."""
    points = dict()
    # look at each package and add PTHs found in library
    for element in board.elements:
        # store element origin
        origin = element.origin
        # alias this footprint
        footprint = board.libraries[element.library][element.footprint]
        # add an entry for each pad in the footprint
        for pad in footprint.smds:
            point = copy.copy(pad.origin)
            if element.rotation != 0:
                point.rotate(element.rotation * math.pi / 180.0)
            if element.mirrored:
                layer = '16'
                point.x = -point.x
            else:
                layer = '1'
            if (element.name, pad.name) in board.contact_refs:
                signal = board.contact_refs[(element.name, pad.name)]
                this_point = (layer, snapped_point(origin + point))
                if signal not in points:
                    points[signal] = set()
                points[signal].add(this_point)
    return points


def read_pad_pths(filename):
    """Read and return pths created by pads from the given file."""
    # hold all pths
    pths = []
    # read footprints from board file
    libraries = read_libraries(filename)
    # read in DRU
    dru = read_design_rules(filename)
    # read element placement in board file
    elements = read_elements(filename)
    # contact references
    contact_refs = read_contact_refs(filename)
    # look at each package and add PTHs found in library
    for element in elements:
        # store element origin
        origin = element.origin
        # alias this footprint
        footprint = libraries[element.library][element.footprint]
        # add an entry for each pad in the footprint
        for pad in footprint.pads:
            point = copy.copy(pad.origin)
            if element.rotation != 0:
                point.rotate(element.rotation * math.pi / 180.0)
            if element.mirrored:
                point.x = -point.x
            outer, inner = pad.get_diameters(dru)
            if (element.name, pad.name) in contact_refs:
                signal = contact_refs[(element.name, pad.name)]
            else:
                signal = None
            pths.append(PTH(origin + point, pad.drill, outer, inner, signal))
    # look at vias and add those
    return pths


def read_design_rules(filename):
    """Return the design result from the given file"""
    # search through the XML tree
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    for dru_xml in root.iter('designrules'):
        return DesignRules(dru_xml)


def read_polygons(filename):
    """Return all polygons in the file."""
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    polygons = []
    for signals in root.iter('signals'):
        for signal in signals.iter('signal'):
            signal_name = signal.attrib['name']
            for polygon in signal.iter('polygon'):
                # skip wires not on layers 1-16
                if int(polygon.attrib['layer']) < 1:
                    continue
                if int(polygon.attrib['layer']) > 16:
                    continue
                polygons.append(Polygon.from_xml(polygon))
                polygons[-1].signal = signal_name
    return polygons


def read_wires(filename):
    """Return all wires in the file."""
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    wires = []
    for signals in root.iter('signals'):
        for signal in signals.iter('signal'):
            signal_name = signal.attrib['name']
            for wire in signal.iter('wire'):
                # skip wires not on layers 1-16
                if int(wire.attrib['layer']) < 1:
                    continue
                if int(wire.attrib['layer']) > 16:
                    continue
                wires.append(Wire.from_xml(wire, signal_name=signal_name))
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
        while True:
            test = (low + high) / 2.0
            point, tangent = wire.get_distance_along(test)
            if test == high or test == low or test == high:
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


def get_nearest_point(point, list_of_points):
    """Return the distance to the nearest point and that point."""
    nearest_distance = None
    nearest_point = None
    for p in list_of_points:
        this_distance = point.distance_to(p)
        if nearest_distance is None or this_distance < nearest_distance:
            nearest_point = p
            nearest_distance = this_distance
    return nearest_distance, nearest_point


def create_teardrop_vias(filename):
    """Create a script to convert pths to teardrop pths in the given file."""
    # hold commands to redraw all pths/wires
    commands = []
    commands.append('set optimizing off;')
    commands.append('change thermals off;')
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
    # hold all pths
    pths = []
    # hold all vias
    vias = read_vias(filename, dru)
    print('Found %d vias' % len(vias))
    # add vias to pths
    for via in vias:
        pths.append(PTH(via.origin,
                        via.drill,
                        via.outer_diameter,
                        via.inner_diameter,
                        via.signal))
    # hold all wires
    wires = read_wires(filename)
    # get wire points
    wire_points = set()
    for wire in wires:
        wire_points.add(wire.p1)
        wire_points.add(wire.p2)
    # read in through hole pads
    if create_teardrops_on_pths:
        pads = read_pad_pths(filename)
        print('Found %d pads' % len(pads))
        pad_tolerance_mm = 10 * native_resolution_mm
        for pad in pads:
            distance, point = get_nearest_point(pad.origin, wire_points)
            if distance > 0 and distance < pad_tolerance_mm:
                pad.origin = point
        pths += pads
        # add through hold pads to pths
    print('Found %d total PTHs' % len(pths))
    # since pad positions are imprecise, snap them to the nearest wire end
    # point if it's close enough
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
    # via points are also fixed
    via_points = set(via.origin for via in pths)
    # store points that are midway through a wire chain
    mid_points = set(key
                     for key, value in wire_count.items()
                     if value == 2 and key[1] not in via_points)
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
    # get map from point to via
    via_at_point = dict()
    for via in pths:
        assert via.origin not in via_at_point
        via_at_point[via.origin] = via
    # hold wire chains which begin at vias
    wire_chains = []
    # hold wires which are used in chains
    used_wires = set()
    for wire in wires:
        layer = wire.layer
        for p in [wire.p1, wire.p2]:
            # point = (wire.layer, p)
            if p not in via_points:
                continue
            used_wires.add(wire)
            # start chain with this wire, set p1 to the via location
            this_chain = [wire]
            wire_chains.append(this_chain)
            if p != wire.p1:
                assert p == wire.p2
                this_chain[0] = wire.reversed()
                assert this_chain[0].p1 == p
            while (layer, this_chain[-1].p2) in mid_points:
                new_wires = wire_by_point[(layer, this_chain[-1].p2)]
                assert len(new_wires) == 2
                if (new_wires[0].p1 == this_chain[-1].p1 or
                        new_wires[0].p2 == this_chain[-1].p1):
                    this_chain.append(new_wires[1])
                    used_wires.add(new_wires[1])
                else:
                    assert (new_wires[1].p1 == this_chain[-1].p1 or
                            new_wires[1].p2 == this_chain[-1].p1)
                    this_chain.append(new_wires[0])
                    used_wires.add(new_wires[0])
                if this_chain[-1].p2 == this_chain[-2].p2:
                    this_chain[-1] = this_chain[-1].reversed()
                else:
                    assert this_chain[-1].p1 == this_chain[-2].p2
    unused_wires = [wire for wire in wires if wire not in used_wires]
    # change wires within each chain such that points are sorted
    # [wire1, wire2, wire3, ...]
    # via at wire1.p1 with wire1.p2 == wire2.p1, etc.
    # for chains which have a via at both ends, delete the second half
    for chain in wire_chains:
        # if there's not a via at the end, keep the entire chain
        if chain[-1].p2 not in via_points:
            continue
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
            point, _ = chain[-1].get_distance_along(alpha)
            chain[-1].p2 = point
            chain[-1].curve *= alpha
    # hold number of teardrops created
    teardrop_count = 0
    # hold wire commands to draw by layer
    wires_by_layer = dict()
    for chain in wire_chains:
        via_point = chain[0].p1
        wire_width = chain[0].width
        via = via_at_point[chain[0].p1]
        via_diameter = via.outer_diameter
        total_chain_length = sum(x.get_length() for x in chain)
        # if chain is too short, don't do anything
        if total_chain_length < via_diameter / 2.0 + teardrop_tolerance_mm:
            continue
        # can't teardrop pths if the wires are bigger than the via diameter
        if via_diameter <= wire_width + teardrop_tolerance_mm:
            continue
        r1 = via_diameter / 2.0
        r2 = teardrop_inner_radius_mm + wire_width / 2.0
        d = math.sqrt((r1 + r2) ** 2 - (r2 + wire_width / 2.0) ** 2)
        result = find_point_on_chain(chain, via_point, d)
        # if chain is not long enough, ignore it
        if result is None:
            # use end of chain
            result = [list(chain[-1].get_distance_along(1.0)),
                      (len(chain) - 1, 1.0)]
            result[0][1] = -result[0][1]
        (junction_point, tangent), (wire_index, alpha) = result
        teardrop_commands = create_teardrop(junction_point,
                                            tangent,
                                            via_point,
                                            (via_diameter - wire_width) / 2.0,
                                            via.signal,
                                            chain[0].width,
                                            create_teardrop_polygons)
        teardrop_count += 1
        # delete portion of chain between via and junction point
        layer = chain[0].layer
        chain[:] = chain[wire_index:]
        if alpha == 1.0:
            chain[:] = chain[1:]
        elif alpha > 0:
            point, _ = chain[0].get_distance_along(alpha)
            chain[0].curve *= 1.0 - alpha
            chain[0].p1 = junction_point
        # for x in teardrop_commands:
        #     print(x)
        if layer not in wires_by_layer:
            wires_by_layer[layer] = []
        wires_by_layer[layer].extend(teardrop_commands)
    # add wires in chains
    for chain in wire_chains:
        if chain[0].layer not in wires_by_layer:
            wires_by_layer[chain[0].layer] = []
        for wire in chain:
            # wires_by_layer[wire.layer].append(wire.get_command(signal_name=wire.signal))
            wires_by_layer[wire.layer].append(wire.get_command())
    # add wires not found in chains
    for wire in unused_wires:
        if wire.layer not in wires_by_layer:
            wires_by_layer[wire.layer] = []
        # wires_by_layer[wire.layer].append(wire.get_command(signal_name=wire.signal))
        wires_by_layer[wire.layer].append(wire.get_command())
    # add wire commands
    for layer in sorted(wires_by_layer.keys()):
        commands.append('layer %s;' % layer)
        commands.extend(sorted(wires_by_layer[layer], key=wire_command_sort))
    print('- Created %d teardrops.' % teardrop_count)
    # set view on top layer
    commands.append('change layer 1;')
    # ratsnest to get rid of airwires
    commands.append('grid last;')
    commands.append('optimize;')
    commands.append('set optimizing on;')
    commands.append('ratsnest;')
    commands.append('group (>0 0);')
    # backup board
    if backup_board_file:
        backup_file(filename)
    # create script filename
    script_filename = os.path.join(os.path.dirname(filename),
                                   'create_teardrops.scr')
    with open(script_filename, 'w') as f:
        f.write('\n'.join(commands))
    print('- Script generated in board file directory.')
    print('- To run, open board and run "script %s".' % script_filename)


def temp():
    board = Board(board_file)
    round_signals(board)
    backup_file(board_file)
    board.generate_script()
    exit(0)
    #create_teardrop_vias(board)


temp()


# execute as a script
if __name__ == "__main__":
    options = sys.argv
    # this_script_filename = options[0]
    del options[0]
    # if no arguments, output docscript
    if not options:
        print(__doc__)
        exit(0)
    # save board filename
    board_filename = options[-1]
    if not os.path.isfile(board_filename):
        print('ERROR: file \"%s\" not found' % board_filename)
    del options[-1]
    # if no options given, output error
    if not options:
        print('ERROR: No options selected.  See usage information.')
        print(__doc__)
        exit(1)
    # get options
    for option in options:
        if option == '--round':
            round_signals(board_filename)
        elif option == '--teardrop':
            create_teardrop_vias(board_filename)
        else:
            print('ERROR: option \"%s\" not recognized' % option)


# snap_wires_to_grid(board_file)

# create_teardrop_vias(board_file)

# delete_via_teardrops(board_file)
