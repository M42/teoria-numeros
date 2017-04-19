#! /bin/python2
# -*- encoding: utf-8 -*-

from sympy import *
from collections import Counter

def libreDeCuadrados(n):
    """Devuelve si el número dado es libre de cuadrados. Espera un entero
    como entrada.
    """
    return all(v == 1 for v in factorint(n).values())


def norma(a, d):
    """Devuelve la norma de un elemento de un cuerpo cuadrático.
    """
    x,y = xy(a,d)
    return simplify((x+y*sqrt(d)) * (x-y*sqrt(d)))

def traza(a, d):
    """Devuelve la traza de un elemento de un cuerpo cuadrático.
    """
    x,y = xy(a,d)
    return simplify(2*x)


def es_entero(a, d):
    """Devuelve si un número es entero algebraico en un cuerpo
    cuadrático. Comprueba que su norma y su traza lo sean.
    """
    return ask(Q.integer(norma(a, d))) and ask(Q.integer(traza(a, d)))


def xy(a,d):
    """Escribe las coordenadas del número del cuerpo algebraico
    en parte racional y coeficiente del radical.
    """
    # Elimina el caso entero, que da error al usar coeff.
    if ask(Q.integer(a)):
        return (a,0)

    # Simplifica el número, necesario para extraer los coeficientes
    a = simplify(a)
    
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
    """Devuelve si alfa es un múltiplo de beta en el anillo cuadrático
    dado por d. Ambos números deben ser enteros algebraicos.
    """

    # Comprueba si el resultado de dividir es un entero
    return es_entero(cociente(a,b,d), d)


def cociente(a,b,d):
    """Devuelve el cociente de la división de alfa por beta en el anillo
    cuadrático dado por d. Ambos números deben ser enteros algebraicos.
    """
    
    # Excepciones en el caso de que no sean enteros algebraicos
    if not es_entero(a,d):
        raise ArithmeticError("{} no es un entero".format(a))
    if not es_entero(b,d):
        raise ArithmeticError("{} no es un entero".format(b))

    # Multiplica por el inverso
    xb,yb = xy(b,d)
    bconjugado = xb - yb*sqrt(d)
    
    return simplify((a * bconjugado) * Rational(1,norma(b,d)))


def eqpell(n,d):
    """Resuelve la ecuación de Pell. Supone d<0.
    """
    # Si n no es cuadrado en módulo d, no puede existir solución a
    # la ecuación x^2 - dy^2 = n.
    if not is_quad_residue(n,-d):
        return []

    # Prueba todos los valores posibles de y, sabiendo que siempre
    # quedará por debajo de la cota y < sqrt(n/-d).
    cota = int(sqrt(n/-d))+1
    solucionespotenciales = ((sqrt(n+d*y**2), y) for y in xrange(cota))
    soluciones = [[(x,y),(x,-y),(-x,y),(-x,-y)] for (x,y) in solucionespotenciales if ask(Q.integer(x))]

    # Devolvemos las soluciones sin repetición
    return list(set([s for sol in soluciones for s in sol]))
        

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
        soluciones = eqpell(4*n,d)
        return [(Rational(x,2)+Rational(y,2)*sqrt(d)) for (x,y) in soluciones if (x-y) % 2 == 0]

    
def es_unidad(a,d):
    """Devuelve si a es una unidad en el anillo de enteros Q(sqrt(d)).
    """
    
    # Retira el caso de que a no sea un entero y ni siquiera
    # esté en el anillo de enteros.
    if not es_entero(a,d):
        return False

    # Un elemento es unidad si y sólo si tiene norma 1 ó -1.
    n = norma(a, d)
    return n == 1 or n == -1


def es_irreducible(a,d):
    """Comprueba si un entero es irreducible.
    """
    # Excepción en el caso de que no sea un entero algebraico
    if not es_entero(a,d):
        raise ArithmeticError("{} no es un entero".format(a))
    
    # Un elemento es irreducible si y sólo si tiene norma número primo
    # o norma número primo al cuadrado y no existe ninǵun entero de
    # norma dicho número primo.
    n = norma(a,d)

    if isprime(n):
        return True
    
    if ask(Q.integer(sqrt(n))) and isprime(sqrt(n)):
        return len(connorma(sqrt(n),d)) == 0

    return False


def factoriza(a,d):
    # Caso base: es una unidad
    if es_unidad(a,d):
        return {}
    
    # Caso base: es irreducible
    if es_irreducible(a,d):
        return {a: 1}
    
    # Prueba a factorizar su norma
    factoresnorma = factorint(norma(a,d))
    factor = factoresnorma.keys()[0]
    
    # Si no existen elementos de norma dada, prueba
    # con el entero primo de esa norma.
    if len(connorma(factor, d)) != 0:
        # Pero si existen elementos de norma dada, prueba con 
        # el primero de ellos que lo divida.
        for factorpotencial in connorma(factor,d):
            factorpotencial = simplify(factorpotencial)
            if divide(a, factorpotencial, d):
                factor = factorpotencial
        
    # Factoriza recursivamente
    dividido = simplify(cociente(a,factor,d))
    return dict(
        Counter({factor: 1}) + 
        Counter(factoriza(dividido, d))
    )
