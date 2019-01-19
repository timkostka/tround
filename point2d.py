import math


class Point2D:
    """Holds a 2D Point class."""

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def dot(self, other):
        return self.x * other.x + self.y * other.y

    def __add__(self, other):
        return Point2D(self.x + other.x, self.y + other.y)

    def __iadd__(self, other):
        self.x += other.x
        self.y += other.y
        return self

    def __sub__(self, other):
        return Point2D(self.x - other.x, self.y - other.y)

    def __isub__(self, other):
        self.x -= other.x
        self.y -= other.y
        return self

    def __mul__(self, other):
        return Point2D(self.x * other, self.y * other)

    def __rmul__(self, other):
        return self * other

    def __imul__(self, other):
        self.x *= other
        self.y *= other
        return self

    def __div__(self, other):
        return Point2D(self.x / other, self.y / other)

    def __itruediv__(self, other):
        self.x /= other
        self.y /= other
        return self

    def __hash__(self):
        return hash((self.x, self.y))

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def normalize(self):
        d = self.norm()
        if d != 0:
            self /= d
        return self

    def __lt__(self, other):
        return self.x < other.x or (self.x == other.x and self.y < other.y)

    def angle(self):
        """Return the angle from the origin to the given point."""
        return math.atan2(self.y, self.x)

    def angle_to(self, other):
        """Return the angle from this point to the given point."""
        return math.atan2(other.y - self.y, other.x - self.x)

    def distance_to(self, other):
        """Return the distance from this point to the given point."""
        return Point2D(self - other).norm()

    def norm(self):
        """Return the distance from the origin."""
        return math.sqrt(self.dot(self))

    def rotate(self, angle):
        """Rotate the point CCW about the origin the given angle."""
        cos = math.cos(angle)
        sin = math.sin(angle)
        x, y = self.x, self.y
        self.x, self.y = x * cos - y * sin, y * cos + x * sin
        return self

    def __repr__(self):
        return '(%g, %g)' % (self.x, self.y)
