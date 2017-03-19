#!/usr/bin/python
# -*- coding: utf-8 -*-

# criptoMRG.py
#
# Implementación de varios algoritmos criptográficos para el curso de
# Teoría de Números y Criptografía de la Universidad de Granada.
#
# Autor: Mario Román García
# Licencia: GPLv3

# Importación de módulos.
#  sympy: Cálculo simbólico.
#  PIL: Manipulación de imágenes.
from sympy import *
from PIL import Image




## PRÁCTICA 1
# Criptografía

## Cifrado RSA
# El cifrado RSA es un cifrado asimétrico usando resultados de teoría
# de números elemental. Se basa en el Teorema de Euler que dice que:
#
#  a^(phi(n)) = 1 mod n
#
# para cualquier a coprimo con n. La asimetría se basa en que calcular
# el phi(n) es computacionalmente difícil cuando no se conoce la
# factorización del número.

def rsa(x,n,e):
    
    """Dada una entrada, módulo y exponente, codifica usando RSA."""
    
    # Usaremos `pow`, que utiliza internamente un algoritmo de
    # exponenciación rápida.
    return pow(x,e,n)

def rsaDec(y,n,phi,e):
    
    """Dado un mensaje codificado en RSA, decodifica conociendo el
    módulo, su totiente y el exponente usados.
    """
    # Calculamos los coeficientes de Bezout
    (u,v,c) = gcdex(e,phi)
    
    # Elevamos para obtener el origen. Nótese que podemos calcular
    # el módulo de u en phi gracias al teorema de Euler.
    return pow(y,long(u)%phi,n)



## Codificación afín de imágenes.
# Codifica afínmente una imagen, con una transformación afín del tipo
# siguiente y tomando módulo 256, para tres canales:
#
#   img(i,j) -> t * img(i,j) + (r,g,b)
#
# Esta es una codificación simétrica. La decodificación es similar y
# simplemente requiere conocer los coeficientes.

def afinimg(image, t, r, g, b):
    
    """Codifica afínmente una imagen dado un coeficiente y las tres
    componentes del vector de desplazamiento.
    """
    
    # Parte la imagen en sus tres canales RGB.
    rchannel, gchannel, bchannel = image.split()

    # Evalúa la codificación afín en sus tres canales
    rnew = Image.eval(rchannel, lambda x: (t*x+r) % 256)
    gnew = Image.eval(gchannel, lambda x: (t*x+g) % 256)
    bnew = Image.eval(bchannel, lambda x: (t*x+b) % 256)

    # Vuelve a unir los tres canales en una sola imagen
    newimage = Image.merge("RGB", (rnew, gnew, bnew))

    return newimage

def afinimgDec(image, t, r, g, b):

    """Decodifica afínmente una imagen dados los coeficientes usados
    para codificar.
    """

    # Calcula la inversa multiplicativa del coeficiente de codificación
    # Comprueba que efectivamente la codificación sea reversible. Si no
    # lo es, lanzará una excepción.
    (tinv, _, gcd) = gcdex(t,256)
    
    if (gcd != 1):
        raise ArithmeticError('The coefficient must be invertible modulo 256!')

    # Codifica utilizando los coeficientes inversos.
    return afinimg(image, tinv, -r*tinv, -g*tinv, -b*tinv)




## Codificación de Vigenère
# La codificación de Vigenère es una codificación simétrica que utiliza
# una matriz como clave.

# Importa el paquete itertools para trabajar con herramientas de la
# programación funcional, que facilitan el tratamiento de listas.
from itertools import cycle

def vigenere(mensaje, matriz, clave):
    
    """Codifica un mensaje con la codificación de Vigenère usando
    la matriz y la clave dadas.
    """
    
    # Busca en la matriz la posición correspondiente al par de letras
    def codificaletra(pair):
        a,b = pair
        c = matriz[0].index(a)
        r = [fila[0] for fila in matriz].index(b)
        return matriz[r][c]
    
    # Codificamos cada letra emparejada con la clave
    cod = map(codificaletra, zip(mensaje,cycle(clave)))
    
    # La presentamos como una cadena en lugar de una lista
    return ''.join(cod) 


def vigenereDec(mensajecodificado, matriz, clave):

    """Decodifica un mensaje con la codificación de Vigenère usando
    la matriz y la clave dadas.
    """

    # Busca en la matriz la posición correspondiente al par de letras
    def decodificaletra(pair):
        c,b = pair
        r = [fila[0] for fila in matriz].index(b)
        d = matriz[r].index(c)
        return matriz[0][d]
    
    # Decodificamos cada letra emparejada con la clave
    decod = map(decodificaletra, zip(mensajecodificado,cycle(clave)))
    
    # La presentamos como cadena en lugar de una lista
    return ''.join(decod)




# PRÁCTICA 2
# Factorización y pseudoprimalidad

def mpot(p, m):
    """Calcula el mayor exponente con el que p divide a m.
    La entrada debe ser un número positivo.
    """
    # Se comenta que en la práctica 10 será necesario que
    # esta función use un número negativo.

    exp = 0
    while (m % p == 0):
        exp = exp+1
        m = m/p

    return exp

def abmod(x,n):
    """Calcula el resto de la división de x por n con signo.  Devuelve el
    representante de menor módulo, que será negativo en la mitad de
    los casos.
    """
    
    return x%n if x < n/2 else x%n - n

def mayorpot(p,x):
    """En el caso p=-1, devuelve 0 si x es no negativo y 1 si es negativo.
    En otro caso, devuelve el exponente de la mayor potencia de p que
    divide a x.

    El valor de p debe ser -1 o positivo.
    """

    if p == -1: return 1 if x < 0 else 0
    return mpot(p, x)

def ssuma(a,b):
    """Devuelve la suma componente a componente de ambas listas.
    En el caso de dos listas de distinta longitud, devuelve None.
    """
    
    # Comprueba que ambas listas tengan la misma longitud
    if len(a) != len(b): return None

    return [x+y for x,y in zip(a,b)]

def parp(l):
    """Cierto si todos los números de la lista son pares.
    """
    return all(x%2 == 0 for x in l)


def bnumer(b, base, n):
    """Comprueba si b es un B-base relativo a n con la
    base dada.
    """
    # Calcula el módulo absoluto del cuadrado
    a = abmod(b**2, n)
    
    # Factoriza el número usando sólo los números de
    # la base. Caso especial en el -1.
    for factor in base:
        if factor != -1:
            a = a / (factor ** mayorpot(factor,a))
        else:
            if a < 0: a = abs(a)
    
    # Devuelve si está completamente factorizado
    return a == 1


def vec_alfa(b, base, n):
    """Comprueba si b es un B-base relativo a n con la
    base dada y en caso de serlo, devuelve su vector alfa.
    """
    
    if not bnumer(b, base, n): return None
    
    a = abmod(b**2,n)
    return [mayorpot(factor,a) for factor in base]


from itertools import combinations
def suma(lista,k):
    return map(lambda l: reduce(lambda x,y: ssuma(x,y),l), combinations(lista,k))
    
def aux(k,r):
    return list(combinations(range(r),k))
