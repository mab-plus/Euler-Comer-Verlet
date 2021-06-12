from modules import deepcopy, np, plt, ticker, cm
from methodes import rayon, Ep, parametres_excentricite, distance_par_revolution

def test_parametres(dt, nb_revolutions, delta, w) :
  # pas de temps
  try :
    dt = float(dt)
    if dt <= 0 or dt > 1 :
      print_usage()
  except ValueError :
    print( '0 < float < 1 en paramètre : %s' % dt, file=sys.stderr )
    print_usage()
  # nb_revolutions
  try :
    nb_revolutions = int(nb_revolutions)
    if nb_revolutions <= 0 :
      print_usage()
  except ValueError :
    print( 'Il faut un entier > = 1 en paramètre : %s' % nb_revolutions, file=sys.stderr )
    print_usage()
  # variation force
  try :
    delta = float(delta)
    if delta == -2. :
      print_usage()
  except ValueError :
    print( 'Il faut un float != -2 en paramètre : %s' % delta, file=sys.stderr )
    print_usage()
  # vent galactique
  try :
    w = float(w)
  except ValueError :
    print( 'Il faut un float en paramètre : %s' % w, file=sys.stderr )
    print_usage()
  return [dt, nb_revolutions, delta, w]

def print_usage() :
  print('usage [arguments] = [méthode dt nb_revolutions]')
  print('       methode                       : Euler, Cromer, Verlet')
  print('       dt (pas de temps)             : 0 < float < 1')
  print('       nb_revolutions              : entier >= 1')
  
  print('usage [arguments] = [méthode dt nb_revolutions k]')
  print('       methode                       : Verlet')
  print('       dt (pas de temps)             : 0 < float < 1')
  print('       nb_revolutions              : entier >= 1')
  print('       k (force du vent galactique)  : float')
  
  print('usage [arguments] = [méthode dt nb_points delta k]')
  print('       methode                       : Verlet, Verlet2')
  print('       dt (pas de temps)             : 0 < float < 1')
  print('       nb_points                   : entier >= 1')
  print('       force en r^(2 + delta)        : float != -2')
  print('       k (force du vent galactique)  : float')
  try :
    exit()
  except :
    exit(-1)


def print_titre(methode) : 
  sep = ''
  s = ''
  for l in range(0, 2*40 + len(methode)) :
    sep = sep + '-'
  for l in range(0, 39) :
    s = s + '-'
  print(sep)
  print(s + ' ' + methode + ' ' + s)
  print(sep)
  return sep

def revolutions(nb_revolutions) :
  if nb_revolutions <= 1 :
    return 'révolution'
  else :
    return 'révolutions'

def print_data_trajectoire(M, nb_points, nb_revolutions) :
  x = M[0]
  y = M[2]
  dM = distance_par_revolution(M, nb_points)
  # Excentricité    
  a, b, excentricite = parametres_excentricite(M, nb_points)

  print('Semi-major axis a, semi-minor axis b and initial eccentricity')
  print('\tM(0)   = ({0}, {1})'.format(x[0], y[0]))
  print('\tM(T/4) = ({0}, {1})'.format(x[nb_points // 4], y[nb_points // 4]))
  print('\ta0     = {0}'.format(a[0]))
  print('\tb0     = {0}'.format(b[0]))
  print('\te0     = {0}'.format(excentricite[0]))
    
  if nb_revolutions > 1 :
    print('Semi-major axis a, semi-minor axis b and eccentricity after {0} '.format(nb_revolutions) + revolutions(nb_revolutions))
    if nb_revolutions == 2 :
      print('\tM(T)        = ({0}, {1})'.format(x[(nb_revolutions-1)*nb_points - 1], y[(nb_revolutions-1)*nb_points - 1])) 
    else :
      print('\tM({0}T)     = ({1}, {2})'.format(nb_revolutions-1, x[(nb_revolutions-1)*nb_points - 1], y[(nb_revolutions-1)*nb_points - 1]))       
    print('\tM({0}T + T/4) = ({1}, {2})'.format(nb_revolutions-1, x[(nb_revolutions-1)*nb_points - 1 + nb_points // 4], y[(nb_revolutions-1)*nb_points - 1 + nb_points // 4])) 
    print('\ta{0}          = {1}'.format(nb_revolutions, a[len(a)-1]))
    print('\tb{0}          = {1}'.format(nb_revolutions, b[len(b)-1]))
    print('\te{0}          = {1}'.format(nb_revolutions, excentricite[len(excentricite)-1]))
    
  if np.mean(excentricite) <= 0.99 :
    if (a[len(a)-1]-b[len(b)-1])/b[len(b)-1] < 1/50 :
      type_trajectoire = 'circular'
    else :
      type_trajectoire = 'elliptical'
  else :
    type_trajectoire = 'opened'
    
  print('\tAverage eccentricity:', np.mean(excentricite))
  print('\tdM = distance(M(0), M({0}T)) = {1:.5}'.format(nb_revolutions, dM[len(dM) - 1]))
  print('orbit seems ' + type_trajectoire)
  return type_trajectoire

def print_data_stabilite_rayon_energie(methode, M, x, y, vx, vy, r, E, dr, dE , nb_points, nb_revolutions, delta) :
  if nb_revolutions > 0 :
    #revolutions pour tester la stabilite de la valeur du rayon et la stabilite de la valeur de l'énergie totale
    # si force en 1/r^2
    print('Stability test of the radius value and the energy value for {0} '.format(nb_revolutions) + revolutions(nb_revolutions) + ' : ')
  else : 
    print('Radius and energy value for {0} points: '.format(nb_points))
  print('\t(x0, y0)   = ({0}, {1})'.format(x[0], y[0]))
  print('\tr0         = {0}'.format(r[0]))
  print('\tdr0        = 0')
  print('\t(vx0, vy0) = ({0}, {1})'.format(vx[0], vy[0]))
  print('\tE0         = {0}'.format(E[0]))
  print(','*(2*40 + len(methode)))  
  
  if nb_revolutions > 0 :
    print('after {0} '.format(nb_revolutions) + revolutions(nb_revolutions))
    N = nb_revolutions
  else :
    print('{0}th point'.format(nb_points))
    N = nb_points
    
  print('\t(x{0}, y{0})   = ({1}, {2})'.format(N, x[len(x)-1], y[len(y)-1]))
  print('\tr{0}         = {1}'.format(N, r[len(r) - 1]))
  print('\t(vx{0}, vy{0}) = ({1}, {2})'.format(N, vx[len(vx)-1], vy[len(vy)-1]))
  print('\tE{0}         = {1}'.format(N, E[len(E) - 1]))
  print(','*(2*40 + len(methode)))
  

  print('Initial, final deviation') 
  print('\tdr = (r{0} - r)/r0    = {1}%'.format(N, dr[len(dr) - 1]))
  print('\tdE = (E{0} - E0)/E0   = {1}%'.format(N, dE[len(dE) - 1]))
  return 1


def parametres_figures_trajectoire(M, r, E, dr, dE, nb_points, nb_revolutions, type_trajectoire, methode, dt, delta, w) :
  date = ''#datetime.datetime.now().strftime('%Y%m%d%H%M%S') 
  dossier = methode + '/'
  x = M[0]
  y = M[2]
  vx = M[1]
  vy = M[3]
  dM = distance_par_revolution(M, nb_points)
  date = ''#datetime.datetime.now().strftime('%Y%m%d%H%M%S') 
  
  if dt < 0.01 :
    delta_t = '{0}'.format(dt)
  else :
    delta_t = '{0:.5e}'.format(dt)
      
  # figure    
  titre = 'Method: ' + methode
  titre = titre + r', $||\overrightarrow{{F}}||\propto \dfrac{{1}}{{  r^{{ {0:1.2f} }} }}$'.format(2 + delta)
  if w == 0 :
    titre = titre + ', no galactic wind'
  else :
    titre = titre + r', galactic wind $\propto{0}GM$'.format(w)

  if methode == 'Verlet_dtau' :
    titre = titre + r', $\delta \tau=$' + delta_t
  else :
    titre = titre + r', $\delta t=$' + delta_t
   
  if nb_revolutions > 0 :
    titre = titre + ', {0} '.format(nb_revolutions) + revolutions(nb_revolutions)
  else :
    if int(nb_points) > 1 :
      titre = titre + ', t={0}'.format(nb_points*dt)
    else :
      titre = titre + ', t={0}'.format(nb_points*dt)
      
  if nb_revolutions == 0 :
    nb_revolutions = nb_points

  label_M0 = r'$M_0$ ({0:.5}, {1:.5})'.format(x[0], y[0]) 
  label_Mf = r'$M_{{  {0}  }}$ ({1:.5}, {2:.5})'.format(nb_revolutions, x[len(x)-1], y[len(y)-1])
  image = dossier + methode + '-' + type_trajectoire + '-' + str(nb_points) + '-' + str(nb_revolutions) + '-' + date + '.png'
  return [titre, label_M0, label_Mf, image]

# figures 
def print_trajectoire(M, d, nb_revolutions, nb_points, titre, label_M0, label_Mf, image) : 
  print(' nb_revolutions',nb_revolutions)
  print(' nb_points',nb_points)
  
  [x,vx,y,vy] = deepcopy(M) 
  
  plt.gca().spines['right'].set_visible(False)
  plt.gca().spines['top'].set_visible(False)
  plt.gca().spines['bottom'].set_position(('data',0))
  plt.gca().spines['left'].set_position(('data',0))
  plt.gca().xaxis.set_ticks_position('bottom')
  plt.gca().yaxis.set_ticks_position('left')
  plt.gca().xaxis.set_major_formatter(plt.NullFormatter())
  plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
  plt.xticks(np.arange(min(x), max(x), 1), color='b', fontsize=8)
  plt.yticks(np.arange(min(y), max(y), 1), color='b', fontsize=8)
  
  if rayon(x[len(x)-1], y[len(y)-1]) > 10*rayon(x[0], y[0]) :
    plt.xlim(-2*x[0], 2*x[0]) 
    plt.ylim(-2*x[0], 2*x[0]) 
  plt.gca().set_aspect("equal", adjustable="box")
    
  # premier point
  plt.plot(x[0], y[0], 'bx', linewidth=0.1, label = label_M0)
  
  #ensembles de couleurs pour les trajectoires
  couleurs = cm.tab20c(np.linspace(0, 1, nb_revolutions))
    
  # trajectoire
  for n in range(1, nb_revolutions + 1) :
    plt.plot(x[nb_points*(n-1) : nb_points*n], y[nb_points*(n-1): nb_points*n], color = couleurs[n-1], lw=1)
  #plt.plot(x, y, 'b', lw=0.2)
  
  #dernier point
  #plt.scatter(x_r[len(x_r)-1], y_r[len(y_r)-1], linewidth=0.5, s=40, facecolors='none', edgecolors='r', label = label_Mf) 
  plt.scatter(x[len(x)-1], y[len(y)-1], linewidth=0.5, s=40, facecolors='none', edgecolors='r', label = label_Mf)
  
  plt.legend(bbox_to_anchor=(1.1, 1.1), loc='best', fontsize=8)
  plt.title(titre, fontsize=8, horizontalalignment='center', y = -0.1)
  #plt.savefig(image, format='png', dpi=200)
  plt.show()
  plt.close()
  return 1


def parametres_figures_stabilite(nb_revolutions, type_trajectoire, methode, dt, delta) :
  date = ''#datetime.datetime.now().strftime('%Y%m%d%H%M%S') 
  dossier = methode + '/'
  if dt < 0.01 :
    delta_t = '{0}'.format(dt)
  else :
    delta_t = '{0:.5e}'.format(dt)
  if delta == 0 :
    titre1 = 'Method ' + methode + ', stability of the radius of an orbit ' + type_trajectoire + r', $\delta t=$' + delta_t + ' pour {0} '.format(nb_revolutions) + revolutions(nb_revolutions) + ', force $\propto \dfrac{{1}}{{  r^{{ {0:1.2f} }} }}$'.format(2+ delta)
    image1 = dossier + 'ecart-rayon-orbite-' + type_trajectoire + '-' + date + '.png'
    titre2 = 'Method ' + methode + ', energy stability of an orbit ' + type_trajectoire + r', $\delta t =$' + delta_t + ' pour {0} '.format(nb_revolutions) + revolutions(nb_revolutions) + ', force $\propto \dfrac{{1}}{{ r^{{ {0:1.2f} }} }}$'.format(2+ delta)
    image2 = dossier + 'ecart-energie-orbite-' + type_trajectoire  + '-' + date + '.png'
    xlabel = 'Revolutions'
  else :
    titre1 = 'Method ' + methode + r', radius variation, $\delta t=$' + delta_t + ', force $\propto \dfrac{{1}}{{  r^{{ {0:1.2f} }} }}$'.format(2+ delta)
    image1 = dossier + 'ecart-rayon-orbite-' + date + '.png'
    titre2 = 'Method ' + methode + r', energy variation, $\delta t =$' + delta_t + ', force $\propto \dfrac{{1}}{{ r^{{ {0:1.2f} }} }}$'.format(2+ delta)
    image2 = dossier + 'ecart-energie-orbite-' + date + '.png'
    xlabel = 'Points'
  return [titre1, image1, titre2, image2, xlabel]


def print_figure_stabilite(ecarts, titre, image, xlabel, ylabel) :
  rv = np.arange(0, len(ecarts), 1)   
  xmin = min(rv)
  xmax = max(rv) + 1
  ymin = min(ecarts)
  ymax = max(ecarts) + 0.1
  
  plt.xlim(xmin, xmax)
  plt.ylim(ymin, ymax)
  plt.plot(np.arange(0, len(ecarts), 1), ecarts, 'r', lw=2)
  
  plt.grid(b=True, which='major', color='b', linestyle=':')
  plt.grid(b=True, which='minor', color='r', linestyle='--', linewidth=0.1)
  
  plt.xlabel(xlabel, fontsize=8)
  plt.ylabel(ylabel, rotation = 'horizontal', fontsize=8)
  plt.title(titre, fontsize=8, horizontalalignment='center', y = 1.05)
  
  plt.tick_params(axis = 'both', labelsize = 6)
 
  plt.gca().spines['top'].set_visible(False)
  plt.gca().spines['right'].set_visible(False)
  plt.gca().xaxis.set_minor_locator(ticker.AutoMinorLocator())
  plt.gca().yaxis.set_minor_locator(ticker.AutoMinorLocator())
  plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter('%1.2f %%'))
  plt.gca().yaxis.set_label_coords(-0.05,1.02)
  
  #plt.savefig(image, format='png', dpi=200)
  plt.show()
  plt.close()
  return 1

