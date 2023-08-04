from typing import List
import numpy as np
import sympy
from scipy.interpolate import lagrange
from numpy.polynomial.polynomial import Polynomial

# Nice https://github.com/j2kun/elliptic-curves-finite-fields
from finitefield.finitefield import FiniteField
from elliptic.elliptic import Ideal
from elliptic.generalized import Point, GeneralizedEllipticCurve
opes = ["+","^","-","*","="]

name_count_mid = 0
variables = []
class Node:
  name = None
  gate = None
  left = None
  right = None

  def __init__(self, name, operand, left, right, node):
    self.name = name
    self.operand = operand
    self.left = left
    self.right = right
    self.node = node

  @classmethod
  def make_node(cls, ope, l, r):
    global name_count_mid
    name_count_mid += 1
    return cls(
      name = "c{}".format(name_count_mid),
      operand = ope,
      left = l,
      right = r,
      node = None,
    )

  @classmethod
  def make_leaf(cls, n):
    if not n.var in variables and not n.is_number():
      variables.append(n.var)
    return cls(
      name = str(n),
      operand = None,
      left = None,
      right = None,
      node = n 
    )

  def is_leaf(self):
    return self.node != None

  def __str__(self):
    if self.node:
      return str(self.node)
    else: 
      return str(self.operand) + "-{}".format(self.name)

class Token:
  num = None
  ope = None
  var = None

  def __init__(self, n, s):
    self.num = n
    if s in opes:
      self.ope = s
    else:
      self.var = s
    
  def is_number(self):
    return self.num != None

  def  __str__(self):
    if self.num:
      return "{}".format(self.num)
    elif self.var: 
      return "{}".format(self.var)
    else: 
      return self.ope

class Gate:
  left = None
  right = None  
  operand = None
  out = None
  def __init__(self, left: str, right: str, operand: str, out: str):
    self.left = left
    self.right = right
    self.operand = operand
    self.out = out 
    
  def __str__(self):
    return "{} {} {} = {}".format(
      self.left,
      self.operand,
      self.right,
      self.out
    )

# x ^ 3 + x + 5 = 35
def tokenize(expr: str) -> List[Token]:
  ret = []
  
  n = ""
  expr += "."
  for i in expr:
    if i == " ":
      pass

    if i.isdigit():
      n += i

    if len(n) != 0:
      if not i.isdigit():
        ret.append(Token(int(n),None))
        ret.append(Token(None,i))
        n = ""
    else:
      ret.append(Token(None,i))

  ret.pop()
  return ret

def parse_pow(tokens: List[Token]):
  if len(tokens) == 1:
    return Node.make_leaf(tokens[0])
  for i in range(len(tokens)):
    if tokens[i].ope == "*":
      left = parse(tokens[0:i])
      right = parse(tokens[i+1:])
      return Node.make_node(tokens[i], left, right)

def parse(tokens: List[Token]):
  if len(tokens) == 1:
    return Node.make_leaf(tokens[0])
  elif len(tokens) == 3:
    l = Node.make_leaf(tokens[0])
    r = Node.make_leaf(tokens[2])
    return Node.make_node(tokens[1], l, r)
  else:
    for i in range(len(tokens)):
      if tokens[i].ope == "+":
        left = parse(tokens[0:i])
        right = parse(tokens[i+1:])
        return Node.make_node(tokens[i], left, right)

    # case like x*x*x
    return parse_pow(tokens)

max_print = 16
gates = []
def print_node(n, dep):
  if not n:
    return

  if type(n) == Token or n.is_leaf():
    print(str(n).rjust(max_print - dep*5))
  else:
    print_node(n.left, dep+1)
    print(str(n).rjust(max_print - dep*5))
    print_node(n.right, dep+1)

def build_gate(n):
  if not n:
    return
  if type(n) == Token:
    return str(n)
  if n.is_leaf():
    return n.name
  else:
    l = build_gate(n.left)
    r = build_gate(n.right)
    gates.append(Gate(l, r, n.operand.ope, n.name))
    return n.name

def calcate(gates, witness):
  var = {}
  for g in gates:
    l = g.left
    if l in witness.keys():
      l = witness[g.left]
    elif l in var:
      l = var[l]
    else:
      l = int(l)
    
    r = g.right
    if r in witness.keys():
      r = witness[g.right]
    elif r in var:
      r = var[r]
    else:
      r = int(r)

    if g.out == "c{}".format(name_count_mid):
      var["out"] = r + l if g.operand == "+" else r * l
    else:
      var[g.out] = r + l if g.operand == "+" else r * l

  return var

def get_v(gate, vector):
  v = [0] * len(vector)
  if gate.left in vector:
    v[vector.index(gate.left)] = 1
  else:
    v[0] = int(gate.left)
  return v

def get_w(gate, vector):
  v = [0] * len(vector)
  if gate.right in vector:
    v[vector.index(gate.right)] = 1
  else:
    v[0] = int(gate.right)
  return v

def get_y(gate, vector):
  v = [0] * len(vector)
  if gate.out in vector:
    v[vector.index(gate.out)] = 1
  else:
    v[-1] = 1
  return v


def caluc_bunbo(t,l):
  ret = 1
  for i in range(1, l+1):
    if i != t:
      ret *= (t - i)
  return ret

def lagrange_polynomial(point, x):
  ret = 0.0
  for j in range(len(point)):
    nume = 1.0
    denom = 1.0
    for i in range(len(point)):
      if i == j: continue
      nume = nume * (x - point[i][0])
      denom = denom * (point[j][0] - point[i][0])
    nume = nume * point[j][1]
    ret = ret + nume / denom
  return ret

def exec(expr, vals, right):
  tokens = tokenize(expr)
  nodes = parse(tokens)
  print_node(nodes, 0)
  print("====")
  print("variables:{}".format(variables))
  vector = [1] + variables + ["c{}".format(x) for x in range(1,name_count_mid)] + ["out"]
  print("vector:{}".format(vector))
  print("vector order: {}".format(len(vector)))
  print("====")
  build_gate(nodes)
  for g in gates:    
    print(g)
  print("====")
  V = np.zeros((0,len(vector)))
  W = np.zeros((0,len(vector)))
  Y = np.zeros((0,len(vector)))
  for gate in gates:
    v = get_v(gate, vector)
    w = get_w(gate, vector)
    y = get_y(gate, vector)
    if gate.operand == "+":
      for i in range(len(vector)):
        v[i] = v[i] + w[i]
      w = [0] * len(vector)
      w[0] = 1
    V = np.append(V, [v], axis=0)
    W = np.append(W, [w], axis=0)
    Y = np.append(Y, [y], axis=0)

  print(V.T)
  print(W.T)
  print(Y.T)

  witness = calcate(gates, vals)
  S = np.array([1] + list(vals.values()) + list(witness.values()))

  print("s", S)
  print("s.V", np.dot(V, S))
  print("s.W", np.dot(W, S))
  print("s.W*s.V",np.dot(V, S)*np.dot(W, S))
  print("s.Y", np.dot(Y, S))
  print("s.W*s.V-s.Y",np.dot(V, S)*np.dot(W, S)-np.dot(Y, S))
  
  Vp = np.zeros((0, len(W.T[0]))) 
  poly_v = []
  for v in V.T:
    x = np.array(range(1,len(v)+1))
    l = lagrange(x, v)
    print("--=")
    print(x)
    print(v)
    print(l)
    print("--")
    poly_v.append(l)
    poly = np.flipud(Polynomial(l.coef[::-1]).coef)
    if len(poly) == 1:
      Vp = np.append(Vp, [np.zeros(len(v))], axis=0)
    else:
      Vp = np.append(Vp, [poly], axis=0)

  Wp = np.zeros((0, len(W.T[0]))) 
  poly_w = []
  for w in W.T:
    x = np.array(range(1,len(w)+1))
    l = lagrange(x, w)
    poly_w.append(l)
    poly = np.flipud(Polynomial(l.coef[::-1]).coef)
    if len(poly) == 1:
      Wp = np.append(Wp, [np.zeros(len(w))], axis=0)
    else:
      Wp = np.append(Wp, [poly], axis=0)

  Yp = np.zeros((0, len(W.T[0]))) 
  poly_y = []
  for y in Y.T:
    x = np.array(range(1,len(y)+1))
    l = lagrange(x, y)
    poly_y.append(l)
    poly = np.flipud(Polynomial(l.coef[::-1]).coef)
    if len(poly) == 1:
      Yp = np.append(Yp, [np.zeros(len(y))], axis=0)
    else:
      Yp = np.append(Yp, [poly], axis=0)

  h_x = np.array(range(1,len(V.T[0])+1))
  h_y = np.zeros(len(W.T[0]))
  print(h_x, h_y)
  H = lagrange(h_x, h_y)
  print()
  print(Vp)
  print(Wp)
  print(Yp)
  print("-@-")
  print(H)

  P = np.dot(Vp.T, S)*np.dot(Wp.T, S)-np.dot(Yp.T, S)
  
  sympy.var('j') 
  t_poly = j - 1 
  for i in range(2,len(P)):
    t_poly *= j - i

  T = sympy.poly(sympy.expand(t_poly)).all_coeffs()
  H = P/T

  # toxic
  t = 10
  k_a = 7
  k_b = 3
  k_c = 5
  k_d = 11

  # https://github.com/amiller/python-zk-proofs/blob/master/simple-zk-proofs.py
  # First the finite field
  q = 2**256 - 2**32 - 2**9 - 2**8 - 2**7 - 2**6 - 2**4 - 1
  Fq = FiniteField(q,1) # elliptic curve over F_q

  # Then the curve, always of the form y^2 = x^3 + {a6}
  curve = GeneralizedEllipticCurve(a6=Fq(7)) # E: y2 = x3+7

  Gx = Fq(0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798)
  Gy = Fq(0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)
  G = Point(curve, Gx, Gy)

  # Proving Key
  # PI_a, PI_a'
  trusted_v = []
  for pv in poly_v:
    trusted_v.append((G*pv(t), G*pv(t)*k_a))

  # PI_b, PI_b'
  trusted_w = []
  for pw in poly_w:
    trusted_w.append((G*pw(t), G*pw(t)*k_b))

  # PI_c, PI_c'
  trusted_y = []
  for py in poly_y:
    trusted_y.append((G*py(t), G*py(t)*k_c))

  # PI_d
  total = []
  for i in range(len(poly_y)):
    total.append(G*(poly_v[i](t) + poly_w[i](t) + poly_y[i](t)) * k_d)

  pi_h = []
  for h in H:
    pi_h.append(G*h)

  gx = []
  for i in range(10):
    gx.append(G*t**i)

  # Vertification Key
  v_ka = G * k_a
  v_kb = G * k_b
  v_kc = G * k_c
  v_z = G * t_poly.subs(j, t)

  return right == witness["out"]

result = exec("x*x*x+x+5", {"x": 3}, 35)
print(result)