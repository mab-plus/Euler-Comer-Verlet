from modules import deepcopy, np,odeint 

# force gravitationnelle
global GM 
GM = 4* np.pi**2

def Fg(r, exp) :
  return -GM/r**exp

# Energie potentiel force Fg
def Ep(r, exp) :
  return GM/((1-exp) * r**(exp - 1))

# vent galactique 
def k(w) :
  return w*GM

def rayon(x, y) :
  return np.sqrt(x**2 + y**2)

# Methodes M=[[x],[vx],[y],[vy]]
def Euler(M, nbre_points, dt) :
  [x,vx,y,vy] = deepcopy(M)
  for n in range(0, nbre_points - 1) :
    #x,y
    x.append(x[n] + vx[n] * dt)
    y.append(y[n] + vy[n] * dt)
    r = rayon(x[n], y[n]) 
    #ax, ay
    ax = Fg(r, 2) * x[n] / r
    ay = Fg(r, 2) * y[n] / r
    #vx, vy
    vx.append(vx[n] + ax * dt)
    vy.append(vy[n] + ay * dt)
  return [x, vx, y, vy]

def Cromer(M, nbre_points, dt) :
  [x,vx,y,vy] = deepcopy(M)
  for n in range(0, nbre_points - 1) :
    r = rayon(x[n], y[n]) 
    #ax, ay
    ax = Fg(r, 2) * x[n] / r
    ay = Fg(r, 2) * y[n] / r
    #vx, vy
    vx.append(vx[n] + ax * dt)
    vy.append(vy[n] + ay * dt)
    #x,y
    x.append(x[n] + vx[n + 1] * dt)
    y.append(y[n] + vy[n + 1] * dt)
  return [x, vx, y, vy]

def Verlet(M, nbre_points, dt, delta, w) :
  [x,vx,y,vy] = deepcopy(M) 
  for n in range(0, nbre_points - 1) :
    r = rayon(x[n], y[n])
    # acceleration n
    ax = Fg(r, 2 + delta) * x[n] / r + k(w)
    ay = Fg(r, 2 + delta) * y[n] / r
    #x,y
    x.append(x[n] + vx[n] * dt + 0.5 * ax * dt**2)
    y.append(y[n] + vy[n] * dt + 0.5 * ay * dt**2)
    r = rayon(x[n+1], y[n+1])
    # acceleration n+1
    ax_p = Fg(r, 2 + delta) * x[n+1] / r + k(w)
    ay_p = Fg(r, 2 + delta) * y[n+1] / r
    #vx, vy
    vx.append(vx[n] + 0.5*(ax + ax_p) * dt)
    vy.append(vy[n] + 0.5*(ay + ay_p) * dt)
  return [x, vx, y, vy]

def Verlet2(M, nbre_points, dt, delta, w) :
  [x,vx,y,vy] = deepcopy(M) 
  for n in range(0, nbre_points - 1) :
    r = rayon(x[n], y[n])
    # acceleration n
    ax = Fg(r, 2 + delta) * x[n] / r + k(w)
    ay = Fg(r, 2 + delta) * y[n] / r
    #x,y
    x.append(x[n] + vx[n] * r**2 * dt + 0.5 * ax * r**4 * dt**2)
    y.append(y[n] + vy[n] * r**2 * dt + 0.5 * ay * r**4 * dt**2)
    r = rayon(x[n+1], y[n+1])
    # acceleration n+1
    ax_p = Fg(r, 2 + delta) * x[n+1] / r + k(w)
    ay_p = Fg(r, 2 + delta) * y[n+1] / r
    #vx, vy
    vx.append(vx[n] + 0.5*(ax + ax_p) * r**2 * dt)
    vy.append(vy[n] + 0.5*(ay + ay_p) * r**2 * dt)
  return [x, vx, y, vy]


def Odeint(M0, nbre_points, dt, delta, w) :
  def f(m, t) :
    x = m[0]
    vx = m[1]
    y = m[2]
    vy = m[3]
    r = rayon(x, y) 
    #ax, ay
    ax = Fg(r, 2 + delta) * x / r + k(w)
    ay = Fg(r, 2 + delta) * y / r
    return [vx, ax, vy, ay]
  
  t = np.linspace(0, dt*nbre_points, nbre_points)
  result = odeint(f, [M0[0][0], M0[1][0], M0[2][0], M0[3][0]], t)
  return [result[:, 0], result[:, 1], result[:, 2], result[:, 3]]


def parametres_orbitaux(M, dt, delta) :
  [x,vx,y,vy] = deepcopy(M)
  # rayon initial et demi grand axe
  r = rayon(x[0], y[0])
  Ec = (vx[0]**2 + vy[0]**2)/2
  E = Ec + Ep(r, 2+delta)
  a = -GM/(2*E)
  # periode suivant la loi de Kepler valable pour une trajectoire elliptique ou circulaire
  #if E < 0 :
  T = np.sqrt(abs(a)**3)
  #else :
    #T = 0
  return T,int(T/dt), E

def parametres_excentricite(M, nbre_points) :
  [x,vx,y,vy] = deepcopy(M)
  a, b, excentricite  =  [ [], [], [] ]
  nbre_revolutions = len(M[0])//nbre_points 
  for k in range(0, nbre_revolutions) :
    ak = max(abs(x[k*nbre_points]), abs(y[k*nbre_points + nbre_points//4]))
    bk = min(abs(x[k*nbre_points]), abs(y[k*nbre_points + nbre_points//4]))
    a.append(ak)
    b.append(bk)
    excentricite.append(np.sqrt(ak**2 - bk**2)/ak) 
  return a, b, excentricite

def distance_par_revolution(M, nbre_points) :
  x = deepcopy(M[0])
  y = deepcopy(M[2])
  dM =  []
  nbre_revolutions = len(x)//nbre_points 
  for k in range(1, nbre_revolutions + 1) :
    dM.append( rayon(x[k*nbre_points-1] - x[0], y[k*nbre_points-1] - y[0]) )
  return dM

def stabilite_rayon_energie(M, nbre_points, nbre_revolutions) :
  [x,vx,y,vy] = deepcopy(M)
  r0 = rayon(x[0], y[0])
  Ec = (vx[0]**2 + vy[0]**2)/2
  E0 = Ec + Ep(r0, 2)
  r, E, dr, dE  = [ [r0], [E0], [0], [0] ]
  #nbre_revolutions = len(M[0])//nbre_points

  for k in range(1, nbre_revolutions + 1) :
    r_k = rayon(x[k*nbre_points-1], y[k*nbre_points-1])
    Ec = (vx[k*nbre_points-1]**2 + vy[k*nbre_points-1]**2)/2
    E_k = Ec + Ep(r_k, 2)
    r.append(r_k)
    E.append(E_k)
    dr.append(100*(r_k - r0) / r0)
    dE.append(abs(100*(E_k -E0)/E0))
  
  return r, E, dr, dE

def stabilite_rayon_energie2(M, nbre_points, delta) :
  [x,vx,y,vy] = deepcopy(M)
  r0 = rayon(x[0], y[0])
  Ec = (vx[0]**2 + vy[0]**2)/2
  E0 = Ec + Ep(r0, 2+delta)
  r, E, dr, dE  = [ [], [], [], [] ]
  #nbre_revolutions = len(M[0])//nbre_points

  for k in range(1, nbre_points) :
    r_k = rayon(x[k-1], y[k-1])
    Ec = (vx[k-1]**2 + vy[k-1]**2)/2
    E_k = Ec + Ep(r_k, 2+delta)
    r.append(r_k)
    E.append(E_k)
    dr.append(100*(r_k - r0) / r0)
    dE.append(abs(100*(E_k -E0)/E0))
  
  return r, E, dr, dE
