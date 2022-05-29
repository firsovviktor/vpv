import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt
from shapely.geometry import Point, Polygon

def integrate(q_1, q_2, p, k, L):
  D = np.abs(q_2 - q_1)
  tau = (q_2 - q_1)/D
  n = (q_2 - q_1)*1j/D
  z = lambda x: q_1 + tau*x - p
  f = lambda x: np.exp(1j * k/(2*L) * np.abs(z(x))**2 ) * (z(x).real*n.real + z(x).imag*n.imag) / (np.abs(z(x))**2) 
  I_real, err = quad(lambda x: np.real(f(x)), 0, D)
  I_imag, err = quad(lambda x: np.imag(f(x)), 0, D)
  #print(I_real, I_imag, ".")
  I = I_real + 1j * I_imag
  return I

def integrate_opening(opening, p, k, L):
  E = 0
  for count in range(opening.vertices.size):
    q_1 = opening.vertices[count]
    next = (count + 1)%opening.vertices.size
    q_2 = opening.vertices[next]
    E += integrate(q_1, q_2, p, k, L)

    if (True):
      eps = 1e-4
      x = np.real(p)
      y = np.imag(p)
      if (np.abs(x)<eps):
        if (np.abs(y)<eps):
          print (integrate(q_1, q_2, p, k, L))
  return E

class Opening:
  pass

openings = []
mode = "polygon"
N = 3
R = 0.001

if (mode == "polygon"):
  ########################
  op1 = Opening()
  #R = 0.001
  op1.vertices = np.array([1+1j]*N)

  phi = 0
  for i in range(N):
    phi -= 2*np.pi/N
    x = R*np.cos(phi)
    y = R*np.sin(phi)
    op1.vertices[i] = x + 1j*y

  vertices_points = []
  for vert in op1.vertices:
    x = np.real(vert)
    y = np.imag(vert)
    p1 = Point(x, y)
    vertices_points.append(p1)
  poly = Polygon(vertices_points)
  op1.polygon = poly
  openings.append(op1)
  #####
  x_min = -1*R
  x_max = R
  x_res = 51

  y_min = -1*R
  y_max = R
  y_res = 51
  ########################
if (mode == "circle"):
  ########################
  op1 = Opening()
  #R = 0.0005
  N = 20
  op1.vertices = np.array([1+1j]*N)

  phi = 0
  for i in range(N):
    phi += 2*np.pi/N
    x = R*np.cos(phi)
    y = R*np.sin(phi)
    op1.vertices[i] = x + 1j*y

  vertices_points = []
  for vert in op1.vertices:
    x = np.real(vert)
    y = np.imag(vert)
    p1 = Point(x, y)
    vertices_points.append(p1)
  poly = Polygon(vertices_points)
  op1.polygon = poly
  openings.append(op1)
  #####
  x_min = -1*R
  x_max = R
  x_res = 60

  y_min = -1*R
  y_max = R
  y_res = 60
  ########################

if (mode == "two-circles"):
  ########################
  op1 = Opening()
  R = 0.5*R
  N = 20
  op1.vertices = np.array([1+1j]*N)

  phi = 0
  for i in range(N):
    phi += 2*np.pi/N
    x = R*np.cos(phi)
    y = R*np.sin(phi)
    op1.vertices[i] = x + 1j*y

  vertices_points = []
  for vert in op1.vertices:
    x = np.real(vert)
    y = np.imag(vert)
    p1 = Point(x, y)
    vertices_points.append(p1)
  poly = Polygon(vertices_points)
  op1.polygon = poly
  openings.append(op1)
  #####
  op1 = Opening()
  R = 20*R
  N = 20
  op1.vertices = np.array([1+1j]*N)

  phi = 0
  for i in range(N):
    phi -= 2*np.pi/N
    x = R*np.cos(phi)
    y = R*np.sin(phi)
    op1.vertices[i] = x + 1j*y

  vertices_points = []
  for vert in op1.vertices:
    x = np.real(vert)
    y = np.imag(vert)
    p1 = Point(x, y)
    vertices_points.append(p1)
  poly = Polygon(vertices_points)
  op1.polygon = poly
  openings.append(op1)
  #####
  x_min = -0.00035
  x_max = 0.00035
  x_res = 60

  y_min = -0.00035
  y_max = 0.00035
  y_res = 60
  ########################

if (False):
  print ("hole")
  for i in range(openings[0].vertices.size):
    print("vert", i, np.real(openings[0].vertices[i]), np.imag(openings[0].vertices[i]))
      

x_points = np.hstack([np.linspace(x_min, x_max, x_res)]*y_res)
y_temp = np.vstack([np.linspace(y_min, y_max, y_res)]*x_res) 
y_points = np.reshape(y_temp.T, -1)
I = np.zeros(x_points.size)

lambd = 7*1e-7
k = 2*np.pi/lambd
L = 0.5
out_counter = -1

for i in range(x_points.size):
  x = x_points[i]
  y = y_points[i]
  p = x + 1j * y
  p1 = Point(x, y)
  E_res = 0
  
  if (True):
    eps = 1e-4
    if (np.abs(x)<eps):
      if (np.abs(y)<eps):
        print (x, y)

  opening_counter = 0
  for opening in openings:
    if (p1.within(opening.polygon)):
      opening_counter += 1 
  if (opening_counter % 2 == 1):
    E_res = 1
  for opening in openings:
    temp_int = integrate_opening(opening, p, k, L)/(2*np.pi)
    E_res -= temp_int
    
    if (True):
      eps = 1e-4
      if (np.abs(x)<eps):
        if (np.abs(y)<eps):
          print ("I", temp_int)
          print ()
          #print (E_res, ".")
          

  out_counter += 1
  if (out_counter%100 == 0):
    print(out_counter)
  I[i] = np.abs(E_res * np.conj(E_res))
  

plt.scatter(x_points, y_points, c=I, cmap = 'viridis')
plt.colorbar()
plt.show()
  
  
  
  
  
