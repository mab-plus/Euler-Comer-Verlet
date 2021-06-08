from modules import *

# force gravitationnelle
global GM 
GM = 4* np.pi**2
#x, vx, y, vy = [1.1, 0., 1., 0.]
#x, vx, y, vy = [1.1, 10., 1., 0.]
#x, vx, y, vy = [1.1, 5., 1., 0.]
x, vx, y, vy = [1.1, 0., 1., -2.]


M0 = [[x], [vx], [y], [vy]]
R1=M0[0]
R2=M0[0]

def print_usage() :
  print('usage [arguments] = [méthode dt nbre_revolutions]')
  print('       methode                       : Euler, Cromer, Verlet, Odeint')
  print('       dt (pas de temps)             : 0 < float < 1')
  print('       nbre_revolutions              : entier >= 1')
  exit()

def energie(vx, vy, r1, r2) :
  #energie cinetiquen potentiel central, potentiel entre planete    #
  Ec = 0.5*(vx**2 + vy**2)
  Ep1 = -GM / r1
  Ep2 = -GM / r2
  return Ec + Ep1 +Ep2

def rayon(x, y) :
  return np.sqrt(x**2 + y**2)

def Euler(M0, nbre_points, dt) :
  [x, vx, y, vy] = deepcopy(M0)
  for n in range(0, nbre_points) :
    #x,y
    x.append(x[n] + vx[n] * dt)
    y.append(y[n] + vy[n] * dt)
    r1 = rayon(x[n], y[n]) 
    r2 = rayon(x[n] - 2, y[n]) 
    #ax, ay
    ax = -GM * ( (x[n] / r1**3) + ( (x[n] - 2) / r2**3) )
    ay = -GM * ( (y[n] / r1**3) + (y[n] / r2**3) )  
    #vx, vy
    vx.append(vx[n] + ax * dt)
    vy.append(vy[n] + ay * dt)
  return [x, vx, y, vy]

def Cromer(M0, nbre_points, dt) :
  [x, vx, y, vy] = deepcopy(M0)
  for n in range(0, nbre_points) :
    r1 = rayon(x[n], y[n]) 
    r2 = rayon(x[n] - 2, y[n]) 
    #ax, ay
    ax = -GM * ( (x[n] / r1**3) + ( (x[n] - 2) / r2**3) )
    ay = -GM * ( (y[n] / r1**3) + (y[n] / r2**3) )  
    #vx, vy
    vx.append(vx[n] + ax * dt)
    vy.append(vy[n] + ay * dt)
    #x,y
    x.append(x[n] + vx[n + 1] * dt)
    y.append(y[n] + vy[n + 1] * dt)
  return [x, vx, y, vy]

def Verlet(M0, nbre_points, dt) :
  [x, vx, y, vy] = deepcopy(M0)
  for n in range(0, nbre_points) :
    r1 = rayon(x[n], y[n]) 
    r2 = rayon(x[n] - 2, y[n]) 
    #ax, ay n
    ax = -GM * ( (x[n] / r1**3) + ( (x[n] - 2) / r2**3) )
    ay = -GM * ( (y[n] / r1**3) + (y[n] / r2**3) )  
    #x,y
    x.append(x[n] + vx[n] * dt + 0.5 * ax * dt**2)
    y.append(y[n] + vy[n] * dt + 0.5 * ay * dt**2)
    r1 = rayon(x[n+1], y[n+1]) 
    r2 = rayon(x[n+1] - 2, y[n+1]) 
    # ax, ay  n+1
    ax_p = -GM * ( (x[n+1] / r1**3) + ( (x[n+1] - 2) / r2**3) )
    ay_p = -GM * ( (y[n+1] / r1**3) + ( y[n+1] / r2**3) ) 
    #vx, vy
    vx.append(vx[n] + 0.5*(ax + ax_p) * dt)
    vy.append(vy[n] + 0.5*(ay + ay_p) * dt)
  return [x, vx, y, vy]


def Odeint(M0, nbre_points, dt) :
  def f(m, t) :
    x, vx, y, vy = deepcopy(m)
    r1 = rayon(x, y) 
    r2 = rayon(x-2, y)
    #ax, ay
    ax = -GM * ( (x / r1**3) + ( (x - 2) / r2**3) )
    ay = -GM * ( (y / r1**3) + (y / r2**3) )  
    return [vx, ax, vy, ay]
  
  t = np.linspace(0, dt*nbre_points, nbre_points)
  result = odeint(f, [M0[0][0], M0[1][0], M0[2][0], M0[3][0]], t)
  return [result[:, 0], result[:, 1], result[:, 2], result[:, 3]]


def ecart_energie_totale(M0, nbre_points) :
  [x, vx, y, vy] = deepcopy(M0)
  r1 = rayon(x[0], y[0]) 
  r2 = rayon(x[0] - 2, y[0]) 
  E0 = energie(vx[0], vy[0], r1, r2)
  ecart = []
  for n in range(0, nbre_points) :
    r1 = rayon(x[n], y[n]) 
    r2 = rayon(x[n] - 2, y[n]) 
    #energie totale
    E = energie(vx[n], vy[n], r1, r2)
    ecart.append(E)
  return ecart, E0, E


def print_figure(x, vx, y, vy, R1, R2, methode) :   
  plt.gca().spines['right'].set_visible(False)
  plt.gca().spines['top'].set_visible(False)
  plt.gca().spines['bottom'].set_position(('data',0))
  plt.gca().spines['left'].set_position(('data',0))
  
  plt.gca().xaxis.set_ticks_position('bottom')
  plt.gca().yaxis.set_ticks_position('left')  
  plt.gca().xaxis.set_major_formatter(plt.NullFormatter())
  plt.gca().yaxis.set_major_formatter(plt.NullFormatter())

  plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(1))
  plt.gca().xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
  plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(1))
  plt.gca().yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
  
  ##R1max = 1.1*max(R1)
  ##R2max =1.1*max(R2)
  ##plt.xlim(-max(R1max, R2max), max(R1max,R2max)) 
  plt.gca().set_aspect("equal", adjustable="box")

  #plt.xticks(fontsize=8)
  #plt.setp(ax.get_xticklabels(), fontsize=8)
 
  # premier point
  plt.scatter(x[0], y[0], s=10, c='purple', marker='x', linewidth=1, edgecolors='purple', clip_on=False, zorder=100, label = 'Planète, départ')
  # trajectoire
  plt.plot(x, y, 'b', lw=0.5, clip_on=False, zorder=0)
  #dernier point
  plt.scatter(x[len(x)-1], y[len(y)-1], s=5, c='purple', linewidth=1, edgecolors='purple', clip_on=False, zorder=100, label = 'Planète, arrivée')

  plt.scatter(0, 0, s=100, c='yellow', linewidth=1, edgecolors='yellow', clip_on=False, zorder=100, label = 'Soleil 1')
  plt.scatter(2, 0, s=100, c='darkorange', linewidth=1, edgecolors='darkorange', clip_on=False, zorder=100, label = 'Soleil 2')

  plt.legend(bbox_to_anchor=(0.9, 1.1), loc='best', fontsize=8)
  plt.title('Method : ' + methode, fontsize=8, horizontalalignment='center', y = -0.1)
  plt.show()
  plt.close()
  
  
def print_ecart_energie_totale(ecart, nbre_points, dt) :   
  t = np.linspace(0, nbre_points*dt, nbre_points)
  xmin = min(t)
  xmax = max(t)
  ymin = min(ecart)
  ymax = max(ecart)
  
  plt.xlim(xmin, xmax)
  plt.ylim(ymin, ymax + 0.02)
  plt.plot(t, ecart, 'r', lw=2)
  
  plt.grid(b=True, which='major', color='b', linestyle=':')
  plt.grid(b=True, which='minor', color='r', linestyle='--', linewidth=0.1)
  
  plt.xlabel('t', fontsize=8)
  plt.ylabel(r'$\frac{{E - E0}}{{E}}\%$', rotation = 'horizontal', fontsize=8)
  plt.title('Method : ' + methode + ', variation d\'énergie', fontsize=8, horizontalalignment='center', y = 1.05)
  
  plt.tick_params(axis = 'both', labelsize = 6)
 
  plt.gca().spines['top'].set_visible(False)
  plt.gca().spines['right'].set_visible(False)
  plt.gca().xaxis.set_minor_locator(ticker.AutoMinorLocator())
  plt.gca().yaxis.set_minor_locator(ticker.AutoMinorLocator())
  plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f %%'))
  
  plt.gca().yaxis.set_label_coords(-0.05,1.02)
  
  plt.show()
  plt.close()
  return 1

########################## main ##########################
if len( sys.argv ) == 1 :
  #print('1',len( sys.argv ))
  print_usage()

# methode
try :
  methode = sys.argv[1]
  if methode != 'Euler' and methode != 'Cromer' and methode != 'Odeint' and methode != 'Verlet' :
    print_usage()
except ValueError :
  print_usage()
  print( '%s' % dt, file=sys.stderr )


# pas de temps
try :
  dt = float(sys.argv[2])
  if dt <= 0 or dt > 1 :
    usage()
except ValueError :
  usage()
  print( '0 < float < 1 en paramètre : %s' % dt, file=sys.stderr )

# nbre_revolutions
try :
   nbre_points = int(sys.argv[3])
   if nbre_points <= 0 :
     usage()
except ValueError :
  usage()
  print( 'Il faut un entier > = 1 en paramètre : %s' % nbre_revolutions, file=sys.stderr )
  
  
print('Initial conditions : ')
print('\tCentral force in 1/r^2') 
print('\tdt         = ', dt)
print('\t(x0, y0)   = ({0}, {1})'.format(x, y))
print('\t(vx0, vy0) = ({0}, {1})'.format(vx, vy))

if methode == 'Euler' :
  x, vx, y, vy = deepcopy(Euler(M0, nbre_points, dt)) 

if methode == 'Cromer' :
  x, vx, y, vy = deepcopy(Cromer(M0, nbre_points, dt))  
  
if methode == 'Verlet' : 
  x, vx, y, vy= deepcopy(Verlet(M0, nbre_points, dt))

if methode == 'Odeint' : 
  x, vx, y, vy= deepcopy(Odeint(M0, nbre_points, dt)) 
  
  
print_figure(x, vx, y, vy, R1, R2, methode) 

M = x, vx, y, vy
ecart, E0, EF = ecart_energie_totale(M,nbre_points)
print_ecart_energie_totale(ecart, nbre_points, dt)

print('\tE0  = {0}'.format(E0))
print('\tE{0} = {1}'.format(int(dt * nbre_points), EF))
