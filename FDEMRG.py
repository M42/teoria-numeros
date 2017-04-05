#! /bin/python2
# -*- encoding: utf-8 -*-

from sympy import *

def libreDeCuadrados(n):
    """Devuelve si el número dado es libre de cuadrados. Espera un entero
    como entrada.
    """
    return all(v == 1 for v in factorint(n).values())


def norma(a):
    """Devuelve la norma de un elemento de un cuerpo cuadrático.
    """
    return a * (a.conjugate())

def traza(a):
    """Devuelve la traza de un elemento de un cuerpo cuadrático.
    """
    return a + (a.conjugate())


def es_entero(a):
    """Devuelve si un número es entero algebraico en un cuerpo
    cuadrático. Comprueba que su norma y su traza lo sean.
    """
    return ask(Q.integer(norma(a))) and ask(Q.integer(traza(a)))


def xy(a,d):
    """Escribe las coordenadas del número del cuerpo algebraico
    en parte racional y coeficiente del radical.
    """
    # Elimina el caso entero, que da error al usar coeff.
    if ask(Q.integer(a)):
        return (a,0)

    # Devuelve los coeficientes
    return (a.coeff(sqrt(d),0), a.coeff(sqrt(d)))

def ab(u,d):
    """Escribe las coordenadas de un entero algebraico en el cuerpo
    algebraico y en la base usual del cuerpo. La entrada debe ser un
    entero algebraico.
    """

    # Obtiene los coeficientes en parte radical y racional y los cambia
    # de base cuando es necesario, esto es, cuando d=1 mod 4.
    x,y = xy(u,d)

    if d % 4 == 1:
        return x-y, y*2
    else:
        return x, y
    
    return (a,b)


def divide(a,b,d):
    """Devuelve si alfa es un divisor de beta en el anillo cuadrático
    dado por d.
    """
    c = (a*b.conjugate())/norma(b,d)
    return es_entero(c)


def cociente(alfa,beta,d):
    """Devuelve el cociente de la división de beta por alfa en el anillo
    cuadrático dado por d.
    """
    pass


def eqpell(n,d):
    """Suponemos d<0.
    """
    # Si n no es cuadrado en módulo d, no puede existir solución a
    # la ecuación x^2 - dy^2 = n.
    if jacobi_symbol(n,-d) == -1:
        return []

    # Prueba todos los valores posibles de y, sabiendo que siempre
    # quedará por debajo de la cota y < sqrt(n/-d).
    cota = int(sqrt(n/-d))+1
    soluciones = ((sqrt(n+d*y**2), y) for y in xrange(cota))
    return [(x,y) for (x,y) in soluciones if ask(Q.integer(sqrt(x)))]
        

def connorma(n,d):
    """Calcula los elementos de O_{d} con norma n.
    """

    if d % 4 != 1:
        # Caso 1: resolvemos la ecuación de Pell
        solucionespell = eqpell(n,d)
        return [x + y*sqrt(d) for (x,y) in solucionespell]
    else:
        # Caso 2: resolvemos la ecucación de Pell
        # teniendo en cuenta la forma de los números
        pass

def es_unidad(a,d):
    """Devuelve si a es una unidad en el anillo de enteros Q(sqrt(d)).
    """
    
    # Retira el caso de que a no sea un entero y ni siquiera
    # esté en el anillo de enteros.
    if not es_entero(a,d):
        return False

    # Un elemento es unidad si y sólo si tiene norma 1 ó -1.
    norma = norma(a)
    return norma == 1 or norma == -1

def es_irreducible(a,d):
    """Comprueba si un número es irreducible.
    """
    # Un elemento es irreducible si y sólo si tiene norma número primo
    # o norma número primo al cuadrado y no existe ninǵun entero de
    # norma dicho número primo.
    norma = norma(a,d)

    if is_prime(norma):
        return True
    
    if is_prime(sqrt(norma)):
        return len(connorma(sqrt(norma),d)) == 0
        
