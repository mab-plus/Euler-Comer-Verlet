from modules import sys, deepcopy, np
from methodes import *
from fn import *
from sys import exit
################################### MAIN ###################################

try :
  data = np.loadtxt('data.txt', delimiter = ';')
except : 
  print('Mauvais fichier de paramètres M0.txt.')
  print_usage()
  
if len(data[0]) != 4 :
  print('Taille de M0 non valable : %s'% M)
  print_usage()
    
# # Conditions initiales, Parametres du fichier data.txt
#[[1.], [0.], [0.], [2*np.pi]    ]  = 1.0;0.0;0.0;6.283185307179586
#[[1.], [0.], [0.], [2*np.pi*1.2]]  = 1.0;0.0;0.0;7.5398223686155035
#[[1.], [0.], [0.], [5.]         ]  = 1.0;0.0;0.0;5.0

# pour choisir les données, 0 pour 1er ligne, 1 2nd ligne etc...
i = 0
M0 = [[data[i][0]], [data[i][1]], [data[i][2]], [data[i][3]]]
x0 = M0[0][0]
y0 = M0[2][0]
vx0 = M0[1][0]
vy0 = M0[3][0]



if len( sys.argv ) == 1 :
  #print('1',len( sys.argv ))
  print_usage()
  
if len( sys.argv ) != 2 and len( sys.argv ) != 4 and len( sys.argv ) != 5 and len( sys.argv ) != 6:
  #print('2',len( sys.argv ))
  print_usage()
  
methode = sys.argv[1]
  
if (methode == 'Euler' or methode == 'Cromer' or methode == 'Verlet') and len( sys.argv ) == 4 :
  dt, nb_revolutions, delta, w = test_parametres(sys.argv[2], sys.argv[3], 0, 0)
elif (methode == 'Verlet'or methode == 'Verlet2') and len( sys.argv ) == 5 :
  dt, nb_revolutions, delta, w = test_parametres(sys.argv[2], sys.argv[3], 0, sys.argv[4])
elif (methode == 'Verlet') and len( sys.argv ) == 6 :
  dt, nb_points, delta, w = test_parametres(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
else :
  print_usage()
  
# ON COMMENCE
sep = print_titre(methode) 

print('Initial conditions: ')
print('\tCentral force in 1/r^{0}'.format(2 + delta)) 

if w != 0 :
  print('\tgalactic wind +{0} x central force'.format(w))
  
print('\t(x0, y0)   = ({0}, {1})'.format(x0, y0))
print('\t(vx0, vy0) = ({0}, {1})'.format(vx0, vy0))
print('\tdt         = ', dt)


# force en 1/r^2
if delta == 0 : 
  T, nb_points, E0 = parametres_orbitaux(M0, dt, delta)
  # on calcul par la méthode demandée
  if methode == 'Euler' :
    M = deepcopy(Euler(M0, nb_points*nb_revolutions, T/nb_points))
  if methode == 'Cromer' :
    M = deepcopy(Cromer(M0, nb_points*nb_revolutions, T/nb_points))
  if methode == 'Verlet' :
    M = deepcopy(Verlet(M0, nb_points*nb_revolutions, T/nb_points, delta, w))
  if methode == 'Verlet2' :
    M = deepcopy(Verlet2(M0, nb_points*nb_revolutions, T/nb_points, delta, w))
    
  # resultats
  x, vx, y, vy = M
  dM = distance_par_revolution(M, nb_points)
  
  print('\tT = {0}'.format(T))
  if int(T/dt) > 1 :
    print('\t{0} points/revolution'.format(int(T/dt)))
  else :
    print('\t{0} point/révolution'.format(int(T/dt)))
  
  print('\tE0 = {0}'.format(E0))
  if E0 < 0 :
    print('The initial total energy is negative, the trajectory must be closed')
  elif E0 == 0 : 
    print('The initial total energy is zero, the trajectory must be open')  
  else  : 
    print('The initial total energy is positive, the trajectory must be open ')

  print(sep)
  type_trajectoire = print_data_trajectoire(M, nb_points, nb_revolutions)
  
  #revolutions pour tester la stabilite de la valeur du rayon et la stabilite de la valeur de l'énergie totale
  print(sep)
  r, E, dr, dE = stabilite_rayon_energie(M, nb_points, nb_revolutions)  
  print_data_stabilite_rayon_energie(methode, M, x, y, vx, vy, r, E, dr, dE , nb_points, nb_revolutions, delta)
  
  titre, labelx, labely, image = parametres_figures_trajectoire(M, r, E, dr, dE, nb_points, nb_revolutions, type_trajectoire, methode, dt, delta, w)
  print_trajectoire(M, delta, nb_revolutions, nb_points, titre, labelx, labely, image)
  
  if nb_revolutions > 1 :
    titre1, image1, titre2, image2, xlabel = parametres_figures_stabilite(nb_revolutions, type_trajectoire, methode, dt, delta)
    print_figure_stabilite(dr, titre1, image1, xlabel, r'$\frac{{r-r_0}}{{r}}$')
    print_figure_stabilite(dE, titre2, image2, xlabel, r'$\left|\frac{{E-E_0}}{{E_0}}\right|$')
  print(sep) 
# force en 1/r^(2+d)
else :
  if int(nb_points) > 1 :
    print('\t{0} points'.format(nb_points))
  else :
    print('\t{0} point'.format(nb_points))
    
  if methode == 'Verlet' :
    M = deepcopy(Verlet(M0, nb_points, dt, delta, w))
  
  # resultats
  x, vx, y, vy = M
  dM = distance_par_revolution(M, nb_points)
  print(sep)

  #revolutions pour tester la stabilite de la valeur du rayon et la stabilite de la valeur de l'énergie totale
  print(sep)
  r, E, dr, dE = stabilite_rayon_energie2(M, nb_points, delta)  
  print_data_stabilite_rayon_energie(methode, M, x, y, vx, vy, r, E, dr, dE , nb_points, 0, delta)
  
  titre, labelx, labely, image = parametres_figures_trajectoire(M, r, E, dr, dE, nb_points, 0, '', methode, dt, delta, w)
  print_trajectoire(M, delta, 1, nb_points, titre, labelx, labely, image)
  
  titre1, image1, titre2, image2, xlabel = parametres_figures_stabilite(25, '', methode, dt, delta)
  print_figure_stabilite(dr, titre1, image1, xlabel, r'$\frac{{r-r_0}}{{r_0}}$')
  print_figure_stabilite(dE, titre2, image2, xlabel, r'$\left|\frac{{E -E_0}}{{E0}}\right|$', )
  print(sep) 


################################### FIN MAIN ###################################

