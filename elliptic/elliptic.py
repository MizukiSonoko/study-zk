
from .generalized import Point

class EllipticCurve(object):
   def __init__(self, a, b):
      # assume we're already in the Weierstrass form
      self.a = a
      self.b = b

      self.discriminant = -16 * (4 * a*a*a + 27 * b * b)
      if not self.isSmooth():
         raise Exception("The curve %s is not smooth!" % self)


   def isSmooth(self):
      return self.discriminant != 0


   def testPoint(self, x, y):
      return y*y == x*x*x + self.a * x + self.b


   def __str__(self):
      return 'y^2 = x^3 + %sx + %s' % (self.a, self.b)


   def __repr__(self):
      return str(self)


   def __eq__(self, other):
      return (self.a, self.b) == (other.a, other.b)


class Ideal(Point):
   def __init__(self, curve):
      self.curve = curve

   def __neg__(self):
      return self

   def __str__(self):
      return "Ideal"

   def __add__(self, Q):
      if self.curve != Q.curve:
         raise Exception("Can't add points on different curves!")
      return Q

   def __mul__(self, n):
      if not (isinstance(n, int) or isinstance(n, long)):
         raise Exception("Can't scale a point by something which isn't an int!")
      else:
         return self

   def __eq__(self, other):
      return type(other) is Ideal

