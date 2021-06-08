from astropy import constants as constantes
from astropy import units as unites

mu = constantes.M_earth * constantes.M_sun /( constantes.M_earth + constantes.M_sun)
print('Masse du soleil                         : %s' % constantes.M_sun)
print('Masse de la terre                       : %s' % constantes.M_earth)
print('Masse reduite                           : %s' % mu)
print('Valeur de la constante gravitationnelle : %s' % constantes.G)
print('Valeur de l\'unité astronomique         : {:.9e} m'.format(unites.AU.to(unites.km)))

