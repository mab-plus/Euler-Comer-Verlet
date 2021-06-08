from modules import *

# force gravitationnelle
global GM 
GM = 4* np.pi**2

x1, vx1, y1, vy1 = [2.52, 0., 0., np.sqrt(GM/2.52)]
x2, vx2, y2, vy2 = [5.24, 0., 0., np.sqrt(GM/5.24)]
M0 = [[x1], [vx1], [y1], [vy1], [x2], [vx2], [y2], [vy2]]


def print_usage() :
  print('usage [arguments] = [méthode dt nbre_revolutions]')
  print('       method                  : Cromer, Odeint')
  print('       dt (step time)          : 0 < float < 1')
  print('       nb_revolutions          : entier >= 1')
  exit()


def energie(vx, vy, r, d, k1, k2) :
  #energie cinetiquen potentiel central, potentiel entre planete    #
  Ec = 0.5*k1*(vx**2 + vy**2)
  Ep1 = -GM * k1 / r
  Ep2 = -GM * k1 * k2 / d
  return Ec + Ep1 +Ep2

def rayon(x, y) :
  return np.sqrt(x**2 + y**2)

def Cromer(M0, nbre_points, dt) :
  [x1, vx1, y1, vy1, x2, vx2, y2, vy2] = deepcopy(M0)
  R1, R2 = [], []
  for n in range(0, nbre_points) :
    r1 = rayon(x1[n], y1[n]) 
    r2 = rayon(x2[n], y2[n]) 
    R1.append(rayon(x1[n], y1[n]) )
    R2.append(rayon(x2[n], y2[n]) )
    r12 = rayon(x2[n] - x1[n], y2[n] - y1[n]) 
    #ax1, ay1
    ax1 = -GM * ( (x1[n] / r1**3) + 0.04*((x2[n] - x1[n])/ r12**3) )
    ay1 = -GM * ( (y1[n] / r1**3) + 0.04*((y2[n] - y1[n])/ r12**3) )  
    #ax2, ay2 
    ax2 = -GM * ( (x2[n] / r2**3) - 0.001*((x2[n] - x1[n])/ r12**3) )
    ay2 = -GM * ( (y2[n] / r2**3) - 0.001*((y2[n] - y1[n])/ r12**3) )
    #vx1, vy1
    vx1.append(vx1[n] + ax1 * dt)
    vy1.append(vy1[n] + ay1 * dt)
    #vx2, vy2
    vx2.append(vx2[n] + ax2 * dt)
    vy2.append(vy2[n] + ay2 * dt)
    #x,y
    x1.append(x1[n] + vx1[n + 1] * dt)
    y1.append(y1[n] + vy1[n + 1] * dt)
    x2.append(x2[n] + vx2[n + 1] * dt)
    y2.append(y2[n] + vy2[n + 1] * dt)
    
  return [x1, vx1, y1, vy1, x2, vx2, y2, vy2, R1, R2]

def Odeint(M0, nbre_points, dt) :
  def f(m, t) :
    x1, vx1, y1, vy1, x2, vx2, y2, vy2 = deepcopy(m)
    r1 = rayon(x1, y1) 
    r2 = rayon(x2, y2)
    r12 = rayon(x2 - x1, y2 - y1)
    #ax1, ay1
    ax1 = -GM * ( (x1 / r1**3) + 0.04*((x2 - x1)/ r12**3) )
    ay1 = -GM * ( (y1 / r1**3) + 0.04*((y2 - y1)/ r12**3) )
    #ax2, ay2 
    ax2 = -GM * ( (x2 / r2**3) - 0.001*((x2 - x1)/ r12**3) )
    ay2 = -GM * ( (y2 / r2**3) - 0.001*((y2 - y1)/ r12**3) )
    return [vx1, ax1, vy1, ay1, vx2, ax2, vy2, ay2]
  
  t = np.linspace(0, dt*nbre_points, nbre_points)
  result = odeint(f, [M0[0][0], M0[1][0], M0[2][0], M0[3][0], M0[4][0], M0[5][0], M0[6][0], M0[7][0]], t)
  return [result[:, 0], result[:, 1], result[:, 2], result[:, 3], result[:, 4], result[:, 5], result[:, 6], result[:, 7]]

def ecart_energie_totale(M, nbre_points) :
  [x1, vx1, y1, vy1, x2, vx2, y2, vy2, r1, r2] = deepcopy(M)
  r12 = rayon(x2[0] - x1[0], y2[0] - y1[0]) 
  E0 = energie(vx1[0], vy1[0], r1[0], r12, 0.001, 0.04) + energie(vx2[0], vy2[0], r2[0], r12, 0.04, 0.001)
  ecart = []
  for n in range(0, nbre_points) :
    r1 = rayon(x1[n], y1[n]) 
    r2 = rayon(x2[n], y2[n]) 
    r12 = rayon(x2[n] - x1[n], y2[n] - y1[n])  
    #energie totale
    E = energie(vx1[n], vy1[n], r1, r12, 0.001, 0.04) + energie(vx2[n], vy2[n], r2, r12, 0.04, 0.001)
    ecart.append(100*(E-E0)/E0)
  return ecart, E0, E


def print_figure(x1, vx1, y1, vy1, x2, vx2, y2, vy2, R1, R2, methode) :     
  plt.gca().spines['right'].set_color('none')
  plt.gca().spines['top'].set_color('none')
  plt.gca().spines['bottom'].set_position(('data',0))
  plt.gca().spines['left'].set_position(('data',0))
  plt.gca().xaxis.set_ticks_position('bottom')
  plt.gca().yaxis.set_ticks_position('left')
  plt.gca().xaxis.set_major_formatter(plt.NullFormatter())
  plt.gca().yaxis.set_major_formatter(plt.NullFormatter())

  R1max = 1.1*max(R1)
  R2max =1.1*max(R2)
  plt.xlim(-max(R1max, R2max), max(R1max,R2max)) 
  plt.gca().set_aspect("equal", adjustable="box")

  # premier point
  plt.plot(x1[0], y1[0], 'bx', linewidth=0.1, label = 'departure planets 1' )
  plt.plot(x2[0], y2[0], 'rx', linewidth=0.1, label = 'departure planets 2' )
  # trajectoire
  plt.plot(x1, y1, 'b', lw=0.2)
  plt.plot(x2, y2, 'r', lw=0.2)
  #dernier point
  plt.scatter(x1[len(x1)-1], y1[len(y1)-1], linewidth=0.5, s=40, facecolors='none', edgecolors='b', label = 'arrivée planètes 1')
  plt.scatter(x2[len(x2)-1], y2[len(y2)-1], linewidth=0.5, s=40, facecolors='none', edgecolors='r', label = 'arrivée planètes 2')

  plt.legend(bbox_to_anchor=(0.9, 1.1), loc='best', fontsize=8)
  plt.title('Méthode : ' + methode, fontsize=8, horizontalalignment='center', y = -0.1)
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
  plt.title('Méthode : ' + methode + ', variation d\'énergie', fontsize=8, horizontalalignment='center', y = 1.05)
  
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
  if methode != 'Cromer' and methode != 'Odeint' :
    print_usage()
except ValueError :
  print_usage()
  print( '%s' % dt, file=sys.stderr )


# pas de temps
try :
  dt = float(sys.argv[2])
  if dt <= 0 or dt > 1 :
    print_usage()
except ValueError :
  print_usage()
  print( '0 < float < 1 in parameter: %s' % dt, file=sys.stderr )

# nbre_revolutions
try :
   nbre_points = int(sys.argv[3])
   if nbre_points <= 0 :
     print_usage()
except ValueError :
  print_usage()
  print( 'You need an integer >= 1 as a parameter: %s' % nbre_revolutions, file=sys.stderr )

  

print('Initial conditions :')
print('\tCentral force in 1/r^2') 
print('\tdt         = ', dt)
print('\t(x1, y1)   = ({0}, {1})'.format(x1, y1))
print('\t(vx1, vy1) = ({0}, {1})'.format(vx1, vy1))
print('\t(x2, y2)   = ({0}, {1})'.format(x2, y2))
print('\t(vx2, vy2) = ({0}, {1})'.format(vx2, vy2))

if methode == 'Cromer' :
  x1, vx1, y1, vy1, x2, vx2, y2, vy2, R1, R2 = deepcopy(Cromer(M0, nbre_points, dt))  

if methode == 'Odeint' : 
  x1, vx1, y1, vy1, x2, vx2, y2, vy2 = deepcopy(Odeint(M0, nbre_points, dt)) 
  R1=M0[0]
  R2=M0[4]
print_figure(x1, vx1, y1, vy1, x2, vx2, y2, vy2, R1, R2, methode) 

M = x1, vx1, y1, vy1, x2, vx2, y2, vy2, R1, R2
ecart, E0, EF = ecart_energie_totale(M,nbre_points)
print_ecart_energie_totale(ecart, nbre_points, dt)

print('\tE0  = {0}'.format(E0))
print('\tE{0} = {1}'.format(int(dt * nbre_points), EF))
