#! /bin/python2
# -*- encoding: utf-8 -*-
from FDEMRG import *
dpositivos = [2,3,5,6,7,11,13,17,19,21,29,33,37,41,57,73]


# es_unidad
# Comprueba distintas unidades.
assert( es_unidad(sqrt(-1),-1) )
assert( es_unidad(1,-1) )

# Comprueba números que no son unidades.
assert( not es_unidad(0,-1) )
assert( not es_unidad(0,3) )
assert( not es_unidad(sqrt(3),3) )


# connorma
# En cuerpos sobre negativos, no podemos encontrar norma negativa.
assert( connorma(-3,-1) == [] )
assert( connorma(-2,-1) == [] )
assert( connorma(-7,-3) == [] )


# pell
# Encuentra una solución de Pell para el caso n=1 para casos conocidos.
for d in dpositivos:
    x,y = pell(d)
    assert( x**2-d*y**2 == 1 )
