#! /bin/python2
# -*- encoding: utf-8 -*-

from random import randint
from sympy import *
from FDEMRG import *

# Enumera los d para los que el anillo cuadrático O(d) es
# un dominio euclídeo.
dpositivos = [2,3,5,6,7,11,13,17,19,21,29,33,37,41,57,73]
dnegativos = [-1,-2,-3,-7,-11]


# libreDeCuadrados
# Comprueba si diferentes números son cuadrados
assert( libreDeCuadrados(7) )
assert( libreDeCuadrados(30) )
assert( libreDeCuadrados(1) )
assert( libreDeCuadrados(-1) )
assert( not libreDeCuadrados(9) )
assert( not libreDeCuadrados(50) )
assert( not libreDeCuadrados(12) )


# norma
# Calcula normas en un anillos cuadrático.
# Normas en cuerpos de enteros negativos.
assert( norma(1+sqrt(-1), -1) == 2 )
assert( norma(1,-1) == 1 )
assert( norma(sqrt(-1), -1) == 1 )
assert( norma(2, -2) == 4 )
assert( norma(1+sqrt(-2), -2) == 3 )
# La norma en anillos cuadráticos negativos
# debería ser siempre positiva
for d in dnegativos:
    a = randint(0,100)
    b = randint(0,100)
    if a + b != 0:
        assert( norma(a + b*sqrt(d), d) > 0 )


# traza
# Calcula traza de elementos en varios cuerpos.
assert( traza(2+sqrt(-2), -2) == 4 )
assert( traza(4-sqrt(-1), -1) == 8 )
assert( traza(3+5*sqrt(3), 3) == 6 )
assert( traza(-1, 5) == -2 )


# es_entero
# Comprueba enteros algebraicos.
assert( es_entero(0,2) )
assert( es_entero(0,-7) )
assert( es_entero(3+sqrt(-1), -1) )
assert( es_entero(2+2*sqrt(41), 41) )


# es_unidad
# Comprueba distintas unidades.
assert( es_unidad(sqrt(-1),-1) )
assert( es_unidad(1,-1) )
assert( es_unidad(3+2*sqrt(2),2) )
assert( es_unidad(1,3) )
# Comprueba números que no son unidades.
assert( not es_unidad(0,-1) )
assert( not es_unidad(0,3) )
assert( not es_unidad(sqrt(3),3) )


# connorma
# En cuerpos sobre negativos, no podemos encontrar norma negativa.
assert( connorma(-3,-1) == [] )
assert( connorma(-2,-1) == [] )
assert( connorma(-7,-3) == [] )
# Encuentra generadores para cada número.
for d in dpositivos:
    generadores = filter(lambda a: a != 1 and a != -1, connorma(1,d))
    assert( generadores != [] )


# pell
# Encuentra una solución de Pell para el caso n=1 para casos conocidos.
for d in dpositivos:
    x,y = pell(d)
    assert( x**2-d*y**2 == 1 )


# eqpell
# Devuelve un generador de la ecuación.
assert( eqpell(1,3) == (2,1) )
assert( eqpell(1,5) == (9,4) )

# idealesDivisores
assert( idealesDivisores(3,3) == [[3,sqrt(3)]] )
assert( idealesDivisores(5,3) == [[5]] )
