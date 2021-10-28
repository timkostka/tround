"""
This script modifies an eagle board file by rounding corners of traces.

Usage:
> python tround.py [options] [board_filename]

For example:
> python tround.py --teardrop /path/to/board.brd
> python tround.py --round /path/to/board.brd


How it works to round traces:

* Read in all wires in the board file.
* Read in all PTHs in the board file.
* Read in all PTHs and SMDs (pads) within components in the board file.
Establish fixed points:
* Rounded wires are not modified and their endpoints are marked as fixed.
* Locations of PTHs are set as fixed
* Locations of SMDs are set as fixed (and axis information recorded).
Find corner and junction points:
* Non-fixed points where 2 wires meet are marked as corners.
* Non-fixed points where 3+ wires meet are marked as junctions.
Find wire chains:
*
Don't round wires near SMD pads:
* We want to avoid the situation where a wire is curved immediately so that it
  intersects other nearby SMDs. Rounding the wire within the SMD pad is not
  desired so, if the wire is long enough, we break the wire in two and fix the
  portion nearest the SMD pad.




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

# if True, will round corners of traces
rounded_corners = True

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
teardrop_inner_radius_mm = 0.080 * 25.4

# if True, will create polygons for teardrops to avoid unplated regions
create_teardrop_polygons = True

# if True, will teardrop vias
create_teardrops_on_vias = True

# if True, will also teardrop plated through holes found in packages
create_teardrops_on_pths = True

# when rounding signals, maximum deviation from the original copper
max_trace_deviation_mils = 10
max_trace_deviation_mils = 10

# when rounding signals, targer inner radius
target_inner_radius_mils = 150
target_inner_radius_mils = 150

# if a corner is rounded less than this distance, ignore it
min_corner_transition_length_mils = 2

# if a junction is rounded less than this distance, ignore it
min_junction_transition_length_mils = 2

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
# board_file = r'\tround\round_traces_test.brd'
# board_file = r'\teardrop_vias\teardrop_test.brd'
# board_file = r'\kct-tester\sandia-cable-tester-rev5-rounded.brd'
# board_file = r'\kct-tester\sandia-cable-tester-rev5.brd'
# board_file = r'\sandia_cable_tester\sandia-cable-tester-rev5.brd'
# board_file = r'\sandia_cable_tester\sandia-cable-tester-rev5-round.brd'
# board_file = r'\tround\round_traces_test_2.brd'
# board_file = r'\sct-adapters\sct-terminal-block-6.brd'
# board_file = r'\sct-adapters\sct-terminal-block-6-round.brd'
# board_file = r'\sct-adapters\sct-bnc-male.brd'
# board_file = r'\micro_ohmmeter\micro_ohmmeter_rev5.brd'
# board_file = r'\octoversity\octoversity_rev3.brd'
board_file = r'\msp430_signal_analysis\test_round_board.brd'

board_file = project_directory + board_file

board_file = r"C:\Users\tdkostk\Documents\msp430_repo\test_board_design\msp430_function_board_rev2.brd"


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

    def get_layer_diameter(self, layer):
        """Return the plated diameter of the PTH on the given layer."""
        if layer == '1' or layer == '16':
            return self.outer_diameter
        assert 1 < int(layer) < 16
        return self.inner_diameter

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
        # self.vertices = [(Point2D(), curve), ...]
        # curve is given in radians
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
                # convert from degrees in file to radians in object
                curve = float(vertex.attrib['curve']) * math.tau / 360.0
            polygon.vertices.append((Point2D(float(vertex.attrib['x']),
                                             float(vertex.attrib['y'])),
                                     curve))
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
            assert -math.tau < curve < math.tau
            command += (' (%s %s) %s' %
                        (format_position(point.x),
                         format_position(point.y),
                         format_angle(curve * 360.0 / math.tau)))
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
        # curve is given in radians
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
        ccw_test = ((ccw_test + math.pi) % math.tau) - math.pi
        cw_test = a1.angle() - theta - a2.angle()
        cw_test = ((cw_test + math.pi) % math.tau) - math.pi
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
            wire.curve = float(wire_xml.attrib['curve']) * math.tau / 360.0
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
        radius = distance / 2.0 / math.sin(self.curve / 2.0)
        radius = abs(radius)
        midpoint = self.p1 + 0.5 * (self.p2 - self.p1)
        normal = (self.p2 - self.p1).rotate(math.tau / 4.0)
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
        angle = start_angle + alpha * self.curve
        point = center + radius * Point2D(math.cos(angle), math.sin(angle))
        tangent = Point2D(math.cos(angle + math.tau / 4.0),
                          math.sin(angle + math.tau / 4.0))
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
        return radius * abs(self.curve)

    def get_angle(self):
        """Return the angle of the wire in degrees."""
        assert self.curve == 0.0
        assert self.get_length() > 0.0
        angle = self.p1.angle_to(self.p2) * 180.0 / math.pi % 360.0
        # snap to 3 decimal places
        angle = round(angle, 3) % 360.0
        return angle

    def get_command(self):
        """Return the command to recreate the wire."""
        if self.curve == 0.0:
            angle = ''
        else:
            assert -math.tau < self.curve < math.tau
            angle = ' %s' % format_angle(self.curve * 360.0 / math.tau)
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

    def generate_script(self):
        """Generate a script which generates copper in this Board."""
        # hold script commands to run
        commands = []
        commands.append('set optimizing off;')
        commands.append('set wire_bend 2;')
        commands.append('grid mm;')
        commands.append('display none;')
        commands.append('display 1 to 16;')
        commands.append('group all;')
        commands.append('ripup (>0 0);')
        commands.append('display none;')
        commands.append('display all;')
        copper_by_layer = dict()
        for wire in self.wires:
            assert type(wire.layer) is str
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
        # create script filename
        script_filename = os.path.join(os.path.dirname(self.filename),
                                       'ttround.scr')
        if verbose:
            print('\nScript generation successful!')
            print('\nTO RUN THE SCRIPT, open the board and run the '
                  'following command:\nscript %s;'
                  % script_filename)
        commands.append('display all;')
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
    """Return the component placements within the file."""
    # tree = ElementTree.parse(filename)
    # root = tree.getroot()
    placements = [Element(element)
                  for element in
                  ElementTree.parse(filename).getroot().iter('element')]
    return placements


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
    # skip very small angles
    if theta * 360.0 / math.tau < 5.0:
        if verbose:
            print('Info: small angle (%g degrees) encountered in '
                  'get_transition_distance() near point %s'
                  % (theta * 360.0 / math.tau, p2))
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
                    layer,
                    polygon_teardrop=False):
    """Return wires and polygons for creating the given teardrop."""
    # hold wires to create for this teardrop
    teardrop_wires = []
    # hold polygons to create for this teardrop
    teardrop_polygons = []
    # tangent should point towards the via
    # (but this is not strictly necessary)
    via = via_point - joining_point
    assert tangent.dot(via) >= 0.0
    # point = point - via_point
    # get direction normal to tangent
    normal = Point2D(tangent.y, -tangent.x)
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
    start_angle_1 = tangent.angle() + math.tau / 4.0
    start_angle_2 = tangent.angle() + math.tau / 4.0
    start_angle_3 = tangent.angle() + math.tau / 4.0
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
    a1 = (end_angle_1 - start_angle_1) * 360.0 / math.tau
    a1 = ((a1 + 180) % 360) - 180
    a2 = (end_angle_2 - start_angle_2) * 360.0 / math.tau
    a2 = ((a2 + 180) % 360) - 180
    a3 = (end_angle_3 - start_angle_3) * 360.0 / math.tau
    a3 = ((a3 + 180) % 360) - 180
    a1 *= math.tau / 360.0
    a2 *= math.tau / 360.0
    a3 *= math.tau / 360.0
    # add initial line joining via to the wire chain
    new_wire = Wire(joining_point, via_point, signal, width, layer, a3)
    teardrop_wires.append(new_wire)
    # add polygon or lines
    if polygon_teardrop:
        new_polygon = Polygon()
        new_polygon.thermals = 'off'
        new_polygon.layer = layer
        new_polygon.rank = '1'
        new_polygon.width = width
        new_polygon.signal = signal,
        new_polygon.vertices.append((joining_point, a1))
        new_polygon.vertices.append((p1, 0.0))
        new_polygon.vertices.append((via_point, 0.0))
        new_polygon.vertices.append((p2, -a2))
        teardrop_polygons.append(new_polygon)
    else:
        new_wire = Wire(joining_point, p1, signal, width, layer, a1)
        teardrop_wires.append(new_wire)
        new_wire = Wire(joining_point, p2, signal, width, layer, a2)
        teardrop_wires.append(new_wire)
    return teardrop_wires, teardrop_polygons


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


def get_smd_length(wire, axes):
    """Return the length required to escape the given SMD axes."""
    if wire.get_length() == 0.0:
        return 0.0
    angle = wire.get_angle() % 180.0
    if angle == axes[2]:
        return axes[0]
    elif (angle + 90) % 180.0 == axes[2]:
        return axes[1]
    return 0.0


def round_signals(board):
    """Round signals within the given Board object."""
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
    # hold points which should not be moved for each signal
    # locked_points[signal_name] = {(layer_name, Point2D, axes), ...}
    locked_points = read_smds_by_signal(board)
    print(locked_points)
    if verbose:
        print('- Found %d SMD pads for %d signals'
              % (sum(len(x) for x in locked_points.values()),
                 len(locked_points)))
    # hold number of corners rounded
    rounded_corner_count = 0
    # hold number of junctions rounded
    rounded_junction_count = 0
    # hold list of all signals to pass through
    signal_names = set(wire.signal for wire in board.wires)
    # sort wires by signal
    wires_by_signal = dict()
    for signal_name in signal_names:
        wires_by_signal[signal_name] = []
    for wire in board.wires:
        wires_by_signal[wire.signal].append(wire)
    # create set of points used by wires in each signal
    signal_wire_points = dict()
    for signal_name, wires in wires_by_signal.items():
        signal_wire_points[signal_name] = set()
        for wire in wires:
            signal_wire_points[signal_name].add(wire.p1)
            signal_wire_points[signal_name].add(wire.p2)
    # TODO: snap SMD points to nearby wires to avoid roundoff issues
    new_locked_points = dict()
    for signal_name in signal_names:
        new_locked_points[signal_name] = set()
        this_locked_points = locked_points[signal_name]
        for value in this_locked_points:
            this_point = value[1]
            if this_point in signal_wire_points:
                new_locked_points[signal_name].add(value)
                continue
            closest_distance = float('inf')
            closest_point = None
            for point in signal_wire_points[signal_name]:
                distance = this_point.distance_to(point)
                if distance < closest_distance:
                    closest_point = point
                    closest_distance = distance
            if closest_distance < 3 * native_resolution_mm:
                print("Snapped SMD point to wire")
                value = (value[0], closest_point, value[2])
            new_locked_points[signal_name].add(value)
    locked_points = new_locked_points
    del new_locked_points
    # hold list of wires to draw to replace all wires corrently on board
    new_wires = []
    # hold list of new polygons
    new_polygons = []
    # loop through each signal
    assert r"N$7" in signal_names
    for signal_name in sorted(signal_names):
        if signal_name == r"N$7":
            print("ASDFASDFASDFASDFA")
            pass
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
        # for wires which have the same start and end point, make it locked
        # for curved wires, make both ends a locked point
        # for wires of different widths, make the corner a locked point
        for wire in these_wires:
            if wire.p1 == wire.p2:
                these_locked_points.add((wire.layer, wire.p1, (0, 0, 0)))
                continue
            if wire.curve != 0.0:
                these_locked_points.add((wire.layer, wire.p1, (0, 0, 0)))
                these_locked_points.add((wire.layer, wire.p2, (0, 0, 0)))
                continue
            for point in [wire.p1, wire.p2]:
                key = (wire.layer, point)
                if key not in corner_width:
                    corner_width[key] = wire.width
                else:
                    if corner_width[key] != wire.width:
                        these_locked_points.add((wire.layer, point, (0, 0, 0)))
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
            #print(these_locked_points)
            for locked_layer, locked_point, _ in these_locked_points:
                #print(locked_layer, locked_point)
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
            if dist < 3 * native_resolution_mm:
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
            if not rounded_corners:
                break
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
            # sort wires by angle from junction to the other point
            wire_angles = []
            for wire in local_wires:
                p2 = wire.get_other_point(point)
                wire_angles.append([(p2 - point).angle(), p2])
            wire_spokes = [(angle, p, wire)
                           for (angle, p), wire in zip(wire_angles,
                                                       local_wires)]
            print(wire_spokes[:5])
            wire_spokes.sort()
            last_wire = list(wire_spokes[0])
            last_wire[0] += math.tau
            wire_spokes.append(tuple(last_wire))
            # get target radius based on our settings
            target_distance = float('inf')
            # target_radius = wire_width + target_inner_radius_mils * 0.0254
            for i in range(len(wire_spokes) - 1):
                p1 = wire_spokes[i][1]
                p3 = wire_spokes[i + 1][1]
                distance = get_transition_distance(wire_width, p1, point, p3)
                # assert distance < 1000
                assert distance > 0.0
                target_distance = min(target_distance, distance)
            # target_distance = target_radius * math.tan(theta / 2.0)
            key = (layer, point)
            assert key not in rounded_distance
            rounded_distance[key] = target_distance
        # adjust SMD information for quicker access
        smd_axes = dict()
        # locked_points[signal_name] = {(layer_name, Point2D, axes), ...}
        for layer, point, axes in locked_points[signal_name]:
            key = (layer, point)
            #assert key not in smd_axes
            smd_axes[key] = axes
        # for wires exiting an SMD pad, adjust rounded distance so it isn't
        # rounded within the SMD pad
        if signal_name == r"N$7":
            print("ASDFASDFASDFASDFA")
        for key, value in rounded_distance.items():
            layer, point = key
            assert key in wires_at_point
            for wire in wires_at_point[key]:
                other_point = wire.get_other_point(point)
                other_key = (layer, other_point)
                if other_key in smd_axes:
                    axes = smd_axes[other_key]
                    wire_length = wire.get_length()
                    allowed_distance = wire_length - get_smd_length(wire, axes)
                    if allowed_distance < 0.0:
                        allowed_distance = 0.0
                    if value > allowed_distance:
                        if super_verbose:
                            print("    - Adjusted rounding distance to "
                                  "avoid SMD pad")
                        rounded_distance[key] = allowed_distance
                        value = allowed_distance
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
            is_corner = len(adjacent_points[key]) == 2
            distance_mils = rounded_distance[key] / 0.0254
            if distance_mils == 0.0:
                del rounded_distance[key]
            elif rounded_distance[key] < native_resolution_mm * 10:
                if super_verbose:
                    print('  - Skipping infinitesimal transition length %g mm'
                          % rounded_distance[key])
                del rounded_distance[key]
            elif (is_corner and
                  distance_mils < min_corner_transition_length_mils):
                if super_verbose:
                    print('  - Skipping corner with small transition length '
                          '%g mm' % rounded_distance[key])
                del rounded_distance[key]
            elif (not is_corner and
                  distance_mils < min_junction_transition_length_mils):
                if super_verbose:
                    print('  - Skipping junction with small transition length '
                          '%g mm' % rounded_distance[key])
                del rounded_distance[key]
        # create remaining part of existing straight segments
        for wire in these_wires:
            # if not modified at all, just copy the wire
            if ((wire.layer, wire.p1) not in rounded_distance and
                    (wire.layer, wire.p2) not in rounded_distance):
                new_wires.append(wire)
                continue
            direction = (wire.p2 - wire.p1).normalize()
            key = (wire.layer, wire.p1)
            pa = copy.copy(wire.p1)
            if key in rounded_distance:
                pa += direction * rounded_distance[key]
            pb = copy.copy(wire.p2)
            key = (wire.layer, wire.p2)
            if key in rounded_distance:
                pb -= direction * rounded_distance[key]
            # if entire wire is rounded, the middle doesn't need drawn
            wire_length = (pa - pb).norm()
            if wire_length < native_resolution_mm:
                continue
            # generate a new wire
            new_wires.append(Wire.from_points(pa,
                                              pb,
                                              wire.width,
                                              wire.signal,
                                              wire.layer))
        # create transitions in corners and junctions
        for (layer, p2), distance in rounded_distance.items():
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
                        angle += math.tau
                    angle = math.pi - angle
                    segments.append((p2a, -angle))
                # create polygon command
                if polygons_in_junctions:
                    polygon = Polygon()
                    polygon.thermals = 'off'
                    polygon.width = corner_width[(layer, p2)]
                    polygon.signal = signal_name
                    polygon.layer = layer
                    polygon.vertices.extend(segments)
                    new_polygons.append(polygon)
                if traces_in_junctions or not polygons_in_junctions:
                    for i in range(len(segments)):
                        new_wire = Wire(segments[i][0],
                                        segments[(i + 1) % len(segments)][0],
                                        signal_name,
                                        corner_width[(layer, p2)],
                                        layer,
                                        segments[i][1])
                        new_wires.append(new_wire)
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
    # add created polygons
    board.new_polygons.extend(new_polygons)


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


def read_pths_by_signal(board, include_vias=True, include_pths=True):
    """Read and return pths by signal for the given board."""
    pths = dict()
    # add vias for each signal
    if include_vias:
        for via in board.vias:
            if via.signal not in pths:
                pths[via.signal] = []
            pths[via.signal].append(PTH.from_via(via))
    # look at each package and add PTHs found in library
    for element in board.elements:
        if not include_pths:
            break
        # store element origin
        origin = element.origin
        # alias this footprint
        footprint = board.libraries[element.library][element.footprint]
        # add an entry for each pad in the footprint
        for pad in footprint.pads:
            point = copy.copy(pad.origin)
            if element.rotation != 0:
                point.rotate(element.rotation * math.tau / 360.0)
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
    """
    Read and return SMD pads by signal for the given board.

    Returned value is a dictionary with a key for each signal name. Each
    dictionary value is a list of the form (layer, point, axes).
    * layer: string value of the layer name
    * point: Point2D value of the center of the SMD pad.
    * axes: triplet of (width, height, angle) where angle is the angle of the
            width axis
    """
    points = dict()
    # look at each package and add PTHs found in library
    for element in board.elements:
        # store element origin
        origin = element.origin
        # alias this footprint
        footprint = board.libraries[element.library][element.footprint]
        # add an entry for each pad in the footprint
        for pad in footprint.smds:
            # set up original axes
            axes = (float(pad.dx) / 2,
                    float(pad.dy) / 2,
                    (pad.rot + element.rotation) % 180.0)
            # get center point
            point = copy.copy(pad.origin)
            if element.rotation != 0:
                point.rotate(element.rotation * math.tau / 360.0)
            if element.mirrored:
                layer = '16'
                point.x = -point.x
            else:
                layer = '1'
            if (element.name, pad.name) in board.contact_refs:
                signal = board.contact_refs[(element.name, pad.name)]
                this_point = (layer, snapped_point(origin + point), axes)
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
                point.rotate(element.rotation * math.tau / 360.0)
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


def teardrop_board_vias(board):
    """Teardrop PTHs and vias within the given Board object."""
    if verbose:
        print('\nCreating teardrops in board')
    # get PTHs by signal
    pths_by_signal = read_pths_by_signal(board,
                                         create_teardrops_on_vias,
                                         create_teardrops_on_pths)
    if verbose:
        print('- Found %d total PTHs in %d signals'
              % (sum(len(x) for x in pths_by_signal), len(pths_by_signal)))
    # get list of all signals within wires and pths
    signal_names = set(wire.signal for wire in board.wires)
    signal_names.update(pths_by_signal.keys())
    # add signal name keys to PTH dict
    for signal in signal_names:
        if signal not in pths_by_signal:
            pths_by_signal[signal] = []
    # sort wires by signal
    wires_by_signal = dict()
    for signal in signal_names:
        wires_by_signal[signal] = []
    for wire in board.wires:
        wires_by_signal[wire.signal].append(wire)
    # hold number of teardrops created
    teardrop_count = 0
    # get list of all wire end points
    wire_points = set()
    for wire in board.wires:
        wire_points.add(wire.p1)
        wire_points.add(wire.p2)
    # hold wires to create
    new_wires = []
    # hold polygons to create
    new_polygons = []
    # since pad positions are imprecise, snap them to the nearest wire end
    # point if it's close enough
    for signal in sorted(signal_names):
        # alias pths for this signal
        these_pths = pths_by_signal[signal]
        # alias wires for this signal
        these_wires = wires_by_signal[signal]
        if super_verbose:
            print('- Processing signal "%s"' % signal)
            print('  - Found %d wires and %d PTHs'
                  % (len(these_wires), len(these_pths)))
        # store wires at each (layer, point)
        wires_at_point = dict()
        for wire in these_wires:
            for point in [wire.p1, wire.p2]:
                this_point = (wire.layer, point)
                if this_point not in wires_at_point:
                    wires_at_point[this_point] = []
                wires_at_point[this_point].append(wire)
        # get map from point to PTH
        pths_at_point = {pth.origin: pth for pth in these_pths}
        # store points that are midway through a wire chain
        mid_points = set(key
                         for key, value in wires_at_point.items()
                         if len(value) == 2 and key[1] not in pths_at_point)
        # store wires used by wire chains
        used_wires = set()
        # store wire chain
        wire_chains = []
        # go through each wire and start wire chains from each PTH
        print('  - Found %d mid points' % len(mid_points))
        #if signal == 'W1':
        #    print(pths_at_point)
        for wire in wires_by_signal[signal]:
            for p in [wire.p1, wire.p2]:
                # point = (wire.layer, p)
                # find closest point
                if not pths_at_point:
                    continue
                if p not in pths_at_point:
                    closest_pth = None
                    closest_distance = float('inf')
                    for point, pth in pths_at_point.items():
                        this_distance = p.distance_to(point)
                        if this_distance < closest_distance:
                            closest_distance = this_distance
                            closest_pth = pth
                    #dist_to_pth = min(
                    #    p.distance_to(x) for x in pths_at_point.keys())
                    if closest_distance > 0.1 * native_resolution_mm:
                        continue
                    print('Found PTH a small distance of %g away' % closest_distance)
                    pths_at_point[p] = closest_pth
                    #if signal == 'W1' and p not in mid_points:
                    #    print('No PTH or midpoint at point %s' % p)
                    #    print('Closest = %g' % dist_to_pth)
                    #continue
                used_wires.add(wire)
                # start chain with this wire, set p1 to the via location
                this_chain = [wire]
                wire_chains.append(this_chain)
                if p != wire.p1:
                    assert p == wire.p2
                    this_chain[0] = wire.reversed()
                    assert this_chain[0].p1 == p
                while (wire.layer, this_chain[-1].p2) in mid_points:
                    this_point = (wire.layer, this_chain[-1].p2)
                    next_wires = wires_at_point[this_point]
                    assert len(next_wires) == 2
                    if (next_wires[0].p1 == this_chain[-1].p1 or
                            next_wires[0].p2 == this_chain[-1].p1):
                        this_chain.append(next_wires[1])
                        used_wires.add(next_wires[1])
                    else:
                        assert (next_wires[1].p1 == this_chain[-1].p1 or
                                next_wires[1].p2 == this_chain[-1].p1)
                        this_chain.append(next_wires[0])
                        used_wires.add(next_wires[0])
                    if this_chain[-1].p2 == this_chain[-2].p2:
                        this_chain[-1] = this_chain[-1].reversed()
                    else:
                        assert this_chain[-1].p1 == this_chain[-2].p2
        # get wires not found in any chain
        unused_wires = [wire for wire in these_wires if wire not in used_wires]
        if super_verbose:
            print('  - Found %d wire chains and %d wires not in chains'
                  % (len(wire_chains), len(unused_wires)))
        # for chains which have a via at both ends, delete the second half
        for chain in wire_chains:
            assert chain[0].p1 in pths_at_point
            # if there's not a via at the end, keep the entire chain
            if chain[-1].p2 not in pths_at_point:
                continue
            # get total length of chain
            chain_length = sum(wire.get_length() for wire in chain)
            assert chain_length > 0.0
            # find approximate midpoint
            start_pth = pths_at_point[chain[0].p1]
            start_length = start_pth.get_layer_diameter(chain[0].layer) / 2.0
            end_pth = pths_at_point[chain[-1].p2]
            end_length = end_pth.get_layer_diameter(chain[-1].layer) / 2.0
            target_length = start_length
            target_length += (chain_length - start_length - end_length) / 2.0
            # find exact midpoint
            length_so_far = 0.0
            index = 0
            while length_so_far + chain[index].get_length() < target_length:
                length_so_far += chain[index].get_length()
                index += 1
            # find proportion along this segment
            alpha = (target_length - length_so_far) / chain[index].get_length()
            if alpha < 0.0:
                alpha = 0.0
            elif alpha > 1.0:
                alpha = 1.0
            # delete trail of chain
            if alpha == 0.0:
                chain[:] = chain[:index]
            else:
                chain[:] = chain[:index + 1]
                if alpha < 1.0:
                    chain[-1] = copy.copy(chain[-1])
                    point, _ = chain[-1].get_distance_along(alpha)
                    chain[-1].p2 = point
                    chain[-1].curve *= alpha
        # create teardrops
        for chain in wire_chains:
            via_point = chain[0].p1
            wire_width = chain[0].width
            pth = pths_at_point[chain[0].p1]
            via_diameter = pth.get_layer_diameter(chain[0].layer)
            total_chain_length = sum(x.get_length() for x in chain)
            # if chain is too short, don't do anything
            if total_chain_length < via_diameter / 2.0 + teardrop_tolerance_mm:
                if super_verbose:
                    print('  - Skipping too short wire chain')
                continue
            # can't teardrop pths if the wires are bigger than the via diameter
            if via_diameter <= wire_width + teardrop_tolerance_mm:
                if super_verbose:
                    print('  - Skipping wire chain with large trace width')
                continue
            r1 = via_diameter / 2.0
            r2 = teardrop_inner_radius_mm + wire_width / 2.0
            d = math.sqrt((r1 + r2) ** 2 - (r2 + wire_width / 2.0) ** 2)
            result = find_point_on_chain(chain, via_point, d)
            # if chain is not long enough, use the last point
            if result is None:
                if super_verbose:
                    print('  - Truncating teardrop because of short trace')
                # use end of chain
                result = [list(chain[-1].get_distance_along(1.0)),
                          (len(chain) - 1, 1.0)]
                distance = via_point.distance_to(result[0][0])
                if distance < (via_diameter + wire_width) / 2.0:
                    if super_verbose:
                        print('  - Truncated chain too short and will be '
                              'skipped')
                    continue
                result[0][1] = -result[0][1]
            (junction_point, tangent), (wire_index, alpha) = result
            radius = (via_diameter - wire_width) / 2.0
            result = create_teardrop(junction_point,
                                     tangent,
                                     via_point,
                                     radius,
                                     chain[0].signal,
                                     chain[0].width,
                                     chain[0].layer,
                                     create_teardrop_polygons)
            teardrop_wires, teardrop_polygons = result
            # if result failed, don't truncate the chain
            if not teardrop_wires and not teardrop_polygons:
                if super_verbose:
                    print('  - Teardrop generation failed')
                continue
            new_wires.extend(teardrop_wires)
            new_polygons.extend(teardrop_polygons)
            teardrop_count += 1
            # delete portion of chain between via and junction point
            chain[:] = chain[wire_index:]
            if alpha == 1.0:
                chain[:] = chain[1:]
            elif alpha > 0:
                point, _ = chain[0].get_distance_along(alpha)
                chain[0].curve *= 1.0 - alpha
                chain[0].p1 = junction_point
        # add wires in chains
        for chain in wire_chains:
            new_wires.extend(chain)
        # add wires not found in chains
        new_wires.extend(unused_wires)
    if verbose:
        print('- Created %d teardrops' % teardrop_count)
    # replace wires with those created
    board.wires = new_wires
    # add new polygons, if any
    board.new_polygons.extend(new_polygons)


def temp():
    board = Board(board_file)
    round_signals(board)
    teardrop_board_vias(board)
    backup_file(board_file)
    board.generate_script()
    exit(0)


# execute as a script
if __name__ == "__main__":
    temp()
    exit(0)
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
            teardrop_board_vias(board_filename)
        else:
            print('ERROR: option \"%s\" not recognized' % option)
