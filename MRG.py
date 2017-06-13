#! /bin/python2
# -*- encoding: utf-8 -*-

# ARCHIVO DE PRÁCTICAS COMPLETAS DE
# TEORÍA ALGEBRAICA DE NÚMEROS Y CRIPTOGRAFÍA.
# Autor: Mario Román García.
# Licencia: GPLv3
#
# Índice de prácticas:
#  1. Criptografía.
#  2. Tests de primalidad.
#  3. Factorización en Dominios Euclideos
#  4. Factorización de ideales
#  5. Número de clase
#

from sympy import *
from PIL import Image
from collections import Counter
from itertools import *
from backports.functools_lru_cache import lru_cache
import numpy as np
import math




##
# PRÁCTICA 1. Criptografía.
##

## Cifrado RSA
# El cifrado RSA es un cifrado asimétrico usando resultados de teoría
# de números elemental. Se basa en el Teorema de Euler que dice que:
#
#  a^(phi(n)) = 1 mod n
#
# para cualquier a coprimo con n. La asimetría se basa en que calcular
# el phi(n) es computacionalmente difícil cuando no se conoce la
# factorización del número. Cuando se conoce la factorización, puede
# calcularse simplemente como
#
#  phi(n) = \prod p_1^{e_1-1}(p_1-i),
#
# donde la factorización de n es \prod p_i^{e_i}.
def rsa(x,n,e):
    """Dada una entrada, módulo y exponente, codifica usando RSA."""
    return pow(x,e,n)

def rsaDec(y,n,phi,e):
    """Dado un mensaje codificado en RSA, decodifica conociendo el módulo,
    su totiente y el exponente usados.
    """
    # Calculamos los coeficientes de Bezout.
    (u,v,c) = gcdex(e,phi)
    # Elevamos para obtener el origen. Nótese que podemos calcular
    # el módulo de u en phi gracias al teorema de Euler.
    return pow(y,long(u)%phi,n)


## Codificación afín de imágenes.
# La codificación afín aplica una transformación de la forma
#  x |-> ax + b
# a todos los datos. Constituye un cifrado simétrico.
#
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
# La codificación de Vigenère es una codificación simétrica que
# utiliza una matriz como clave. Dada esa matriz, simplemente busca la
# fila y columna determinadas por la letra a codificar y la clave que
# le corresponde y usa la letra en esa entrada como
# sustituta. Normalmente se repite la clave tantas veces como sea
# necesario, como en
#
#  MENSAJECIFRADO
#  CLAVECLAVECLAV,
#
# la inversa se calcula trivialmente.

def vigenere(mensaje, matriz, clave):    
    """Codifica un mensaje con la codificación de Vigenère usando la
    matriz y la clave dadas.
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
    """Decodifica un mensaje con la codificación de Vigenère usando la
    matriz y la clave dadas.
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













##
# PRÁCTICA 2. Factorización y pseudoprimalidad
##
# Primalidad: un elemento a es primo si a|bc implica a|b o a|c.
# Irreducibilidad: un elemento a es irreducible si b|a implica b~a ó b~1.
#
# En todo DFU coinciden primalidad e irreducibilidad.

def mpot(p, m):
    """Calcula el mayor exponente con el que p divide a m.
    La entrada debe ser un número positivo.
    """
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
    return x%n if x%n < n/2 else x%n - n

def mayorpot(p,x):
    """En el caso p=-1, devuelve 0 si x es no negativo y 1 si es negativo.
    En otro caso, devuelve el exponente de la mayor potencia de p que
    divide a x.

    El valor de p debe ser -1 o positivo.
    """
    if p == -1: return 1 if x < 0 else 0
    return mpot(p, x)

def ssuma(a,b):
    """Devuelve la suma componente a componente de ambas listas. En el
    caso de dos listas de distinta longitud, devuelve None.
    """
    # Comprueba que ambas listas tengan la misma longitud
    if len(a) != len(b): return None

    return [x+y for x,y in zip(a,b)]

def parp(l):
    """Cierto si todos los números de la lista son pares.
    """
    
    return all(x%2 == 0 for x in l)


def bnumer(b, base, n):
    """Comprueba si b es un B-base relativo a n con la base dada.
    """
    # Calcula el módulo absoluto del cuadrado
    a = abmod(b**2, n)
    
    # Factoriza el número usando sólo los números de
    # la base. Caso especial en el -1.
    for factor in base:
        a = a / (factor ** mayorpot(factor,a))
    
    # Devuelve si está completamente factorizado
    return a == 1


def vec_alfa(b, base, n):
    """Comprueba si b es un B-base relativo a n con la base dada y en caso
    de serlo, devuelve su vector alfa.
    """
    if not bnumer(b, base, n): return None
    
    a = abmod(b**2,n)
    return [mayorpot(factor,a) for factor in base]

def suma(lista,k):
    return map(lambda l: reduce(lambda x,y: ssuma(x,y),l), combinations(lista,k))
    
def aux(k,r):
    return list(combinations(range(r),k))

def primerosNPrimos(n):
    """Devuelve una lista con los primeros n números primos.
    """
    return [prime(i) for i in xrange(1,n)]


## Test 1: Primalidad
# Pseudoprimalidad
# Comprueba que n es pseudoprimo en base b
def isPsp(n,b):
    return gcd(b,n) == 1 and pow(b,n-1,n) == 1

# Test simple utilizando una sola base aleatoria
def pspTest(n):    
    # Elige una base aleatoria
    b = random.randint(2,n-1)
    if not isPsp(n,b):
        return b, False
    return b, True

# Test de pseudoprimalidad completo
def psp(n,k=1):
    # Prueba el test de primalidad simple sobre un conjunto de k bases
    # elegidas aleatoriamente.
    bases = []
    for _ in xrange(k):
        # Si alguna falla, sabemos inmediatamente que no es primo.
        b,result = pspTest(n)
        if result == False: return (b,False)
        bases.append(b)
    
    # Si todas las bases pasan el test, es posible que sea primo.
    print "Es posible que {} sea primo".format(n)
    return bases,True


## Test 2: Pseudoprimalidad de Euler
# Símbolo de Jacobi.
# El símbolo de Jacobi se define como:
# 
# jacobi(a,n) = \prod legendre(a,p_i)^{\alpha_i},
# 
# donde los p_i dan la factorización de n y donde legendre representa
# el símbolo de Legendre, siendo legendre(a,p)
#
#  * 0 cuando a = 0 mod p.
#  * 1 cuando a es residuo cuadrático mod p.
#  * -1 cuando a no es residuo cuadrático mod p.
#
# El cálculo del símbolo de Jacobi lo realizamos sabiendo las siguientes
# propiedades
#
#  * jacobi(m,n) * jacobi(n,m) = (-1)^{(m-1)/2 * (n-1)/2}
#  * jacobi(m,n) = jacobi(m+kn,n)
#  * jacobi(2,n) = (-1)^{(n^2-1)/8}
#  * jacobi(1,n) = 1
#
# y usando recursión.
def jacobiSym(m,n):
    """Calcula el símbolo de Jacobi de dos enteros. El segundo de ellos
    no debe ser par."""
    
    # Trabaja en módulo n
    m = m % n
    
    # Caso base, no son coprimos
    if gcd(m,n) != 1:
        return 0
    
    # Caso base, 1 es un residuo cuadrático 
    if m == 1:
        return 1
    
    # Caso recursivo, m divisible entre 2
    if m % 2 == 0:
        return pow(-1, (n**2-1)/8) * jacobiSym(m/2, n)
    
    # Caso recursivo, ley de reciprocidad cuadrática
    return pow(-1, ((m-1)*(n-1))/4) * jacobiSym(n,m)

# Pseudoprimalidad de Euler
# Definición de pseudoprimos de Euler.
def isEpsp(n,b):
    """Comprueba la pseudoprimalidad de Euler de un
    número en una base dada.
    """
    
    # Comprueba coprimalidad
    if gcd(b,n) != 1:
        print "{} es divisor de {}".format(gcd(b,n), n)
        return False
    
    # Test de Euler
    return pow(b,(n-1)/2,n) == jacobi_symbol(b,n)


# Test simple de pseudoprimalidad con una base aleatoria
# sólo para números impares.
def epspSimpleTest(n):
    """Prueba la pseudoprimalidad de Euler de un número en
    una base aleatoria. Devuelve el número y si ha funcionado
    el test.
    """
    
    # Elige una base aleatoria
    b = random.randint(2,n-1)
    if not isEpsp(n,b):
        return b, False

    return b, True

# Test de pseudoprimalidad de Euler
def epsp(n,k=1):
    """Prueba el test de pseudoprimalidad de Euler con k bases
    distintas aleatorias.
    """
    # Eliminamos el caso de un número par
    if n % 2 == 0:
        print "{} es par".format(n)
        return 2,False
    
    # Prueba el test de primalidad simple sobre un conjunto de k bases
    # elegidas aleatoriamente.
    bases = []
    for _ in xrange(k):
        
        # Si alguna falla, sabemos inmediatamente que no es primo.
        b, isepsp = epspSimpleTest(n)
        if not isepsp: return (b,False)
        bases.append(b)
    
    # Si todas las bases pasan el test, es posible que sea primo.
    print "Es posible que {} sea primo".format(n)
    return bases,True


## Test 3: fuertemente pseudoprimos.
# Un número n es fuertemente pseudoprimo respecto de b si
#
#  * n es impar,
#  * gcd(b,n) = 1, y
#  * n = 2^st nos da b^t=1, o b^{t2^i} = -1 mod n.
#
# Test para fuertemente pseudoprimos.
def fpsp(n, k=1):
    
    # Comprueba el test de primalidad para una base y una división
    # de n en partes par e impar dadas
    def fpspTest(n,b,s,t):
        # Comprueba coprimalidad
        if gcd(b,n) != 1:
            print "{} es divisor de {}".format(gcd(b,n),n)
            return b, False

        # Comprueba la primera condición de pseudoprimalidad de Euler
        powb = pow(b,t,n)
        if powb == 1 or powb == -1: return b, True

        # Comprueba la segunda condición de pseudoprimalidad de Euler
        for i in xrange(1,s):
            if pow(b, t*2**i, n) == n-1:
                return b, True
            if pow(b, t*2**i, n) == 1:
                return b, False
        
        # Si falla la comprobación
        return b, False
    
    # Caso unidad
    if n == 1: return 1, False
    
    # Caso par
    if n % 2 == 0:
        print "{} es par".format(n)
        return 2,n == 2
    
    # Descompone al número en parte par e impar
    s = mpot(2, n-1)
    t = (n-1)/s
    
    # Elige una base al azar para cada paso
    bases = []
    for _ in xrange(k):
        b = randint(2,n-1)
        bases.append(b)
        b, test = fpspTest(n,b,s,t)
        if not test: return b, False
    
    # En caso de que haya pasado todos los tests
    return bases, True

## Factor-Base
# Elige una lista de B-números con la esperanza de que sean suficientes
# para resolver x^2 = y^2 mod n.
def bi(n,k,i,base):
    """Elige una lista de B-números dados índices sobre los que buscar
    y una base.
    """
    # Crea las listas
    lista1 = [int(math.floor(math.sqrt(i*n))) for i in xrange(1,k+1)]
    lista2 = [m+i for i in xrange(i) for m in lista1]
    
    # Selecciona los B-números
    return filter(lambda b: bnumer(b, base, n), lista2)

# Resolución de la ecuación x^2 = y^2 mod n.
def soleqBase(n, base, bnumeros, k, i):
    # Calcula alfa vectores de los B-números dados.
    alfavectores = map(lambda b: vec_alfa(b,base,n), bnumeros)
    
    # Busca en el conjunto potencia de los B-números. Para cada subconjunto,
    # calcula sus alfa-vectores y comprueba si suman un exponente par.
    # Para las pares, intenta encontrar una solución no trivial.
    powerset = chain.from_iterable(
        combinations(range(len(alfavectores)), a)
        for a in range(len(alfavectores)+1)
    )
    
    for tupla in powerset:
        if tupla != ():
            suma = reduce(lambda x,y: ssuma(x,y),map(lambda x: alfavectores[x], tupla))
            if parp(suma):
                # Calcula la solución, comprueba que no es trivial
                t = reduce(lambda x,y: (x*y)%n, map(lambda x: bnumeros[x], tupla))
                s = reduce(lambda x,y: (x*y)%n, map(lambda (p,e): pow(p,(e/2),n), zip(base, suma)))
                if t != s and t != (-s)%n:
                    return t%n,s%n

def soleq(n, h, k, i):
    # Escoge una base y crea B-números en esa base.
    base = [-1] + primerosNPrimos(h)
    bnumeros = bi(n,k,i,base)
    return soleqBase(n, base, bnumeros, k, i)

## Factorización usando la resolución de Factor-Base.
def fac(n,h,k=5,i=5):
    factorization = {}
    
    # Caso unidad
    if n == 1:
        return {}
    
    # Caso de potencias de 2 y números pares
    if mpot(2,n) != 0:
        factorization[2] = mpot(2,n)
        return factorization
    
    # Caso primo
    if isprime(int(n)):
        factorization[n] = 1
        return factorization
    
    # Caso compuesto
    a,b = soleq(int(n),h,k,i)
    u,v = gcd(a+b, n), n/gcd(a+b, n)
    factorization = dict(Counter(fac(u,h,k,i))+Counter(fac(v,h,k,i)))
    
    return factorization

def soleqfc(n,s):
    # Cálculo de las fracciones continuas y numeradores
    f = continued_fraction_periodic(0,1,n)
    l1 = [f[0]]+f[1] if len(f)>1 else [f[0]]
    l2 = continued_fraction_convergents(l1[:s])
    pbs = [fraction(x)[0] for x in l2]
    
    # Factorización de los cuadrados en módulo absoluto
    factorizaciones = [factorint(abmod(b**2,n)) for b in pbs]
    
    # Base factor
    # Primera condición, aparecen para al menos dos b's
    from collections import Counter
    counter = Counter(flatten(map(lambda d: d.keys(), factorizaciones)))
    baseFactor = [k for k,v in counter.iteritems() if v > 1]
    counter2 = flatten([[k for k,v in d.iteritems() if v%2 == 0] for d in factorizaciones])
    counter2 = [k for k,v in Counter(counter2).items() if v==1]
    baseFactor = baseFactor + counter2
    baseFactor
    
    # Calcula los B-números
    bumerosObtenidos = filter(lambda x: bnumer(x,baseFactor,n), pbs)
    
    return soleqBase(n, baseFactor, bumerosObtenidos, 1, 1)

def facfc(n, s=145):
    factorization = {}
    
    if n == 1:
        return {}
    
    # Caso cuadrado perfecto
    # Para este caso no encuentra solución el algoritmo anterior
    sqrtn = int(sqrt(n))
    if sqrtn != 1 and n % sqrtn == 0:
        return dict(Counter(facfc(sqrtn))+Counter(facfc(sqrtn)))
                    
    # Caso de potencias de 2 y números pares
    if mpot(2,n) != 0:
        factorization[2] = mpot(2,n)
        return factorization
    
    # Caso primo
    if isaPrime(int(n)):
        factorization[n] = 1
        return factorization
    
    # Caso compuesto
    a,b = soleqfc(int(n),s)
    u,v = gcd(a+b, n), n/gcd(a+b, n)
        
    factorization = dict(Counter(facfc(u))+Counter(facfc(v)))
    
    return factorization






















##
# PRÁCTICA 3. Factorización en dominios euclídeos
##
# Llamamos cuerpos de números a una extensión de Q que sea finita.
# Toda extensión de Q finita tiene un elemento primitivo; en particular
# toda extensión cuadrática es de la forma Q(sqrt(d)) para algún d libre
# de cuadrados.
#
# En estos cuerpos de números, llamamos enteros algebraicos a aquellos
# números que son raíz de un mónico con coeficientes en Z.
#
#  1. a es entero algebraico ssi Irr(a) tiene coeficientes en Z.
#  2. a entero algebraico si N(a),T(a) son enteros.
#  3. un racional es entero algebraico ssi es entero.
#  4. sumas y productos de enteros algebraicos son enteros algebraicos.
#  5. los enteros algebraicos forman un anillo.
#  6. para cualquier elemento a existe un n tal que na es entero algebraico.
#  7. todo cuerpo de números está generado con un entero algebraico.
#  8. los enteros algebraicos son un grupo abeliano finitamente generado.
#
# El anillo de enteros algebraicos es Dominio Euclídeo para
#   d = -1,-2,-3,-7,-11,
# y es un Dominio de Factorización Única para
#   d = -19,-43,-67,-163.

def libreDeCuadrados(n):
    """Devuelve si el número entero dado es libre de cuadrados. Espera un
    entero como entrada, devuelve un booleano.

    """
    # Será libre de cuadrados si todos los factores de su
    # factorización tienen exponente 1.
    return all(v == 1 for v in factorint(n).values())


def norma(a, d):
    """Devuelve la norma de un elemento de un cuerpo cuadrático. Toma como
    entrada el elemento del cuerpo cuadrático y la d que define el
    cuerpo cuadrático O(d). Devuelve un natural en la salida para los
    enteros algebraicos.
    """
    x,y = xy(a,d)
    return simplify((x + y*sqrt(d)) * (x - y*sqrt(d)))


def traza(a, d):
    """Devuelve la traza de un elemento de un cuerpo cuadrático.  La traza
    se define como la suma del elemento con su conjugado, es decir,
    como el doble de su parte racional.
    """
    x,y = xy(a,d)
    return simplify(2*x)


def es_entero(a, d):
    """Devuelve si un número es entero algebraico en un cuerpo
    cuadrático. Comprueba que su norma y su traza lo sean.
    """
    return sympify(norma(a, d)).is_integer and sympify(traza(a, d)).is_integer


def xy(a,d):
    """Escribe las coordenadas del número del cuerpo algebraico en parte
    racional y coeficiente del radical.
    """
    # Elimina el caso entero, que da error al usar coeff.
    a = expand(a)
    if sympify(a).is_integer:
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
    """Devuelve si alfa es un múltiplo de beta en el anillo cuadrático
    dado por d. Ambos números deben ser enteros algebraicos.
    """
    # Comprueba si el resultado de dividir es un entero
    return es_entero(cociente(a,b,d), d)


def cociente(a,b,d):
    """Devuelve el cociente de la división de alfa por beta en el anillo
    cuadrático dado por d. Ambos números deben ser enteros
    algebraicos.

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


def pell(d):
    """Resuelve la ecuación de Pell para n=1 y d>0.
    """
    # Calcula la fracción continua y el numerador y denominador
    # de la última fracción de los convergentes.
    f = continued_fraction_periodic(0,1,d)
    l = flatten(f)[:-1]
    convergentes = continued_fraction_convergents(l)
    ultimo = list(convergentes)[-1]
    x0, y0 = ultimo.as_numer_denom()
    
    # Separa según la paridad de la lista de coeficientes
    if (len(l) % 2 == 0):
        return x0, y0
    else:
        return x0*x0+d*y0*y0, 2*x0*y0


def generalpell(n,d):
    """Resuelve la ecuación de Pell para d<0.
    """
    # Resuelve primero para el caso n=1.
    r,s = pell(d)
    if n == 1:
        return [(r,s),(r,-s),(-r,s),(-r,-s)]
    
    # Marca las cotas donde buscará las soluciones.
    cotainf = 0
    cotasup = 0
    if n > 0:
        cotainf = 0
        cotasup = sqrt(Rational(n*(r-1), 2*d))
    elif n < 0:
        cotainf = sqrt(Rational(-n,d))
        cotasup = sqrt(Rational(-n*(r+1), 2*d))
    
    # Encuentra las soluciones de la ecuación de Pell.
    solucionespotenciales = ((n+d*y**2, y) for y in xrange(cotainf,cotasup+1))
    soluciones = (
        [(x,y),(x,-y),(-x,y),(-x,-y)]
        for (x,y) in ((sqrt(x),y) for (x,y) in solucionespotenciales if esCuadrado(x))
    )
    
    return (s for sol in soluciones for s in sol)


cuadrados = set([0,1,4,9,16,25,36,49,64,81,100,121,144,169,196,
 225,256,289,324,361,400,441,484,529,576,625,676,
 729,784,841,900,961,1024,1089,1156,1225,1296,1369,
 1444,1521,1600,1681,1764,1849,1936,2025,2116,2209,
 2304,2401,2500])

def esCuadrado(n):
    """Devuelve si un entero es un cuadrado perfecto."""
    # Los únicos cuadrados en módulo 16 terminan en 0,1,4,9. Así,
    # podemos filtrar unos cuantos números al empezar.
    # n = int(n)
    hx = n & 0xF
    if (hx > 9): return False
    if (hx!=2 and hx!=3 and hx!=5 and hx!=6 and hx!=7 and hx!=8):
        if n in cuadrados:
            return True
        if n < 2500:
            return False
        
        t = math.floor(math.sqrt(n) + 0.5)
        return t*t == n
        
    return False

def eqpell_neg(n,d):
    """Resuelve la ecuación de Pell para d<0.
    """
    n = int(n)
    
    # Si n es negativo, como d también lo es, no puede existir
    # solución a la ecuación
    if n < 0:
        return []
    
    # Si n no es cuadrado en módulo d, no puede existir solución a
    # la ecuación x^2 - dy^2 = n.
    if not is_quad_residue(n,-d):
        return []
    
    # Prueba todos los valores posibles de y, sabiendo que siempre
    # quedará por debajo de la cota y < sqrt(n/-d).
    cota = int(sqrt(n/-d))+1
    solucionespotenciales = ((n+d*y**2, y) for y in xrange(cota))
    soluciones = (
        [(x,y),(x,-y),(-x,y),(-x,-y)]
        for (x,y) in ( (sqrt(x),y) for (x,y) in solucionespotenciales if esCuadrado(x) )
    )

    # Devolvemos las soluciones sin repetición
    return list(set([s for sol in soluciones for s in sol]))


def eqpell(n,d):
    """Resuelve la ecuación de Pell. En el caso de que d sea positivo sólo
    devolverá las soluciones generadoras.
    """
    n = int(n)
    if n == 1 and d > 0:
        return pell(d)
    
    if d <= 0:
        return eqpell_neg(n,d)
    else:
        return generalpell(n,d)

@lru_cache(maxsize=5)
def connorma(n,d):
    """Calcula los elementos de O_{d} con norma n.
    """
    if d % 4 != 1:
        # Caso 1: resolvemos la ecuación de Pell
        solucionespell = eqpell(n,d)
        return (x + y*sqrt(d) for (x,y) in solucionespell)
    
    else:
        # Caso 2: resolvemos la ecucación de Pell
        # teniendo en cuenta la forma de los números
        soluciones = eqpell(4*n,d)
        return (Rational(x,2)+Rational(y,2)*sqrt(d) for (x,y) in soluciones if (x-y) % 2 == 0)
    
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
    
    # Comprueba si la raíz de un número es primo. Si lo es, intenta
    # comprobar si ramifica buscando soluciones.
    raiz = sqrt(n)
    if ask(Q.integer(raiz)) and isprime(raiz):
        return len(list(connorma(raiz,d)) + list(connorma(-raiz,d))) == 0
    
    return False


def factoriza(a,d):
    """Factoriza el elemento en el cuerpo O(d).
    """
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
    connormadada = chain(connorma(factor, d), connorma(-factor,d))
    # Pero si existen elementos de norma dada, prueba con 
    # el primero de ellos que lo divida.
    for factorpotencial in connormadada:
        factorpotencial = simplify(factorpotencial)
        if divide(a, factorpotencial, d):
            factor = factorpotencial
            break
        
    # Factoriza recursivamente
    dividido = simplify(cociente(a,factor,d))
    return dict(
        Counter({factor: 1}) + 
        Counter(factoriza(dividido, d))
    )


def e(d):
    """Devuelve el e del anillo de enteros algebraicos. Según el caso
    será uno de
    
      * 1/2 + 1/2 sqrt(d), cuando d%4 == 1.
      * sqrt(d),           cuando d%4 != 1.
    """
    if d % 4 == 1:
        return Rational(1,2)+Rational(1,2)*sqrt(d) 
    else:
        return sqrt(d)


def irr_e(d):
    """Devuelve el polinomio irreducible de e en los racionales.
    """
    if d % 4 == 1:
        def irr(x):
            return x**2 - x + (1-d)/4
        return irr
    else:
        def irr(x):
            return x**2 - d
        return irr


def irreducible_e(d):
    """
    Devuelve el polinomio irreducible de e en la extensión
    dada por O(d). Este polinomio irreducible tiene como
    coeficientes la traza y la norma.
    """
    return [1, -traza(e(d)), norma(e(d))]









##
# PRÁCTICA 4. Factorización de ideales
##

def matrizrel(ideal, d):
    """Calcula la matriz de relatores asociada a un ideal. La matriz de
    relatores tiene en sus columnas la lista de generadores del ideal
    como grupo abeliano. Con ella tendremos una presentación del
    grupo, que puede simplificarse por transformaciones elementales.
    
    El ideal se introduce como su lista de generadores como ideal.  La
    matriz de relatores se devuelve como lista de listas.
    """
    # Calcula los generadores del ideal como grupo abeliano y para cada
    # generador, calcula sus coordenadas en la base.
    generadores = [x for par in map(lambda a: [a, expand(e(d)*a)], ideal) for x in par]
    return map(lambda a: list(ab(a,d)), generadores)


def LR(matriz):
    """Calcula la forma reducida asociada a una matriz. La forma reducida
    será la matriz equivalente triangular inferior.
    """
    # Ordena la matriz por valor absoluto, dejando los ceros al final.
    matriz = sorted(matriz,
                    key=lambda x: abs(x[0]) if x[0] != 0 else float("inf"))
    
    # Si ha simplificado ya la primera fila, simplifica sólo la
    # segunda fila.
    if (len(matriz) < 2 or matriz[1][0] == 0):
        matriz = [matriz[0]] + sorted(matriz[1:],
                                      key=lambda x: abs(x[1]) if x[1] != 0 else float("inf"))
        
        if (len(matriz) < 3 or matriz[2][1] == 0):
            # Si ya ha simplificado segunda fila, ha terminado
            return matriz
        
        a = matriz[1][1]
        for i in range(2,len(matriz)):
            matriz[i][1] = matriz[i][1] % a
        
        return LR(matriz)
        
    # Si no la ha simplificado, reduce desde el pivote.
    a = matriz[0][0]
    b = matriz[0][1]
    for i in range(1,len(matriz)):
        coc = floor(matriz[i][0] / a)
        matriz[i][0] = matriz[i][0] - a * coc
        matriz[i][1] = matriz[i][1] - b * coc
    
    return LR(matriz)


def simplificaIdeal(ideal, d):
    """Simplifica un ideal para escribirlo como generado por, a lo sumo,
    dos generadores.
    """
    # Toma la matriz de relatores reducida
    reducida = LR(matrizrel(ideal,d))
    
    # Caso particular de longitud 1 o 0. El ideal está ya simplificado
    if len(reducida) < 2:
        return reducida

    # Toma los dos primeros elementos de la matriz de relatores y
    # los escribe como elementos del ideal.
    return [reducida[0][0] + reducida[0][1]*e(d), reducida[1][1]*e(d)]


def normaIdeal(ideal, d):
    """Calcula la norma de un ideal. Viene dado como lista de
    generadores.
    """
    # Reduce la matriz de relatores del ideal. La norma se obtiene
    # como el producto de la primera subdiagonal.
    reducida = LR(matrizrel(ideal, d))
    return abs(reducida[0][0] * reducida[1][1])


def esO(ideal, d):
    """Comprueba si un ideal es el total. Esto es, si es igual a O. El
    ideal viene dado como lista de generadores.
    """
    # Un ideal es de enteros sólo si todos sus generadores lo son.
    # Un ideal de enteros es el total si y sólo si tiene norma 1.
    return all([es_entero(e,d) for e in ideal]) and normaIdeal(ideal, d) == 1 


def pertenece(u, ideal, d):
    """Comprueba si el elemento pertenece o no a un ideal dado.
    """
    # Intentaremos resolver el sistema diofántico:
    #     ax      = n
    #     bx + cy = m
    # donde (n,m) son las coordenadas del elemento y (a,b,c) los
    # elementos de la matriz reducida.
    
    # Primero creamos la matriz reducida de relatores del ideal. De
    # ella tomamos los elementos necesarios para plantear el sistema
    # de ecuaciones diofánticas del que buscaremos solución. Tomamos
    # también las coordenadas del elemento.
    reducida = LR(matrizrel(ideal, d))
    a = reducida[0][0]
    b = reducida[0][1]
    c = reducida[1][1]
    n,m = ab(u,d)
    
    # Resolvemos x si es posible
    if n % a == 0:
        x = n/a
    else:
        return False
    
    # Resolvemos y si es posible
    m = m - x * b
    return m % c == 0


def divideIdeal(I, J, d):
    """Comprueba si el ideal J divide al ideal I. Ambos se introducen como
    listas de generadores.
    """
    # Un ideal divide a otro si y sólo si lo contiene. Y un ideal
    # contiene a otro si contiene a todos sus generadores.
    return all([pertenece(i,J,d) for i in I])


def idealesDivisores(p,d):
    """Dado un número primo, busca los ideales irreducibles que lo
    dividen; es decir, que dividen al ideal generado por ese número
    primo.
    """
    # Analizamos si el primo ramifica o no considerando el polinomio
    # irreducible de e en el cuerpo de característica el primo dado.
    # Ramificará si y sólo si existen raíces a ese polinomio.
    irr = irr_e(d)
    solucionesirr = [x for x in xrange(p) if irr(x)%p==0]
    ramifica = (len(solucionesirr) != 0)
    
    if ramifica:
        # Cuando ramifica, tenemos divisores dados por las soluciones
        # al polinomio característico.
        return [ [p,e(d)-a] for a in solucionesirr]
    else:
        # Para los primos que no ramifican, sabemos que pO será primo
        # y con norma p^2.
        divisor = [p]
        return [divisor]


def productodos(I,J,d):
    """Calcula el producto de dos ideales. Ambos se introducen como
    listas de generadores. Devuelve una lista con dos generadores.
    """
    # Puede calcularse el producto de dos ideales uniéndolos y
    # obteniendo su matriz reducida. Los generadores serán los que
    # queden en esa matriz reducida.
    return simplificaIdeal([expand(i*j) for i in I for j in J],d)


def producto(listaideales,d):
    """Calcula el producto de una lista de ideales. Devuelve dos
    generadores del ideal producto.
    """
    return reduce(lambda x,y: productodos(x,y,d), listaideales, [1,e(d)])


def factoriza_id(ideal, d):
    """Factoriza un ideal. La salida del algoritmo será una lista de
    ideales que serán primos o unidades.
    
        * Si es el total, hemos terminado la descomposición.
        * Para comprobar si es primo, calculamos la norma del ideal y
            * si su norma es primo, es un ideal primo.
            * si es el cuadrado de un primo, será primo sólo si el primo
              no ramifica. Tenemos que estudiar si O/(p)O es un cuerpo o no.
              Para ello calculamos el irreducible de e; p ramifica si y sólo
              si f_p, el polinomio a módulo p, es irreducible.
        * Si la norma del ideal no es primo, la factorizamos y tomamos uno
          de sus primos para calcular sus divisores primos. Si no ramifica, 
          su ideal generado es válido.
        * Si ramifica, se descompone como dos ideales de norma p_i, que serán
          precisamente <p_1,e-a> y <p_1,e-b> (se da sin demostración) para los
          a,b en los que ha factorizado el polinomio irreducible en módulo p
          como:  f_p = (X-a)(X-b).
        * Para cada divisor encontrado, lo añadimos a la lista de factorización
          y seguimos trabajando con el cociente.
    """
    # Factoriza sólo ideales
    #assert( type(ideal) == list )
    # No puede factorizar el ideal cero
    #assert( ideal != [] )
    #assert( normaIdeal(ideal,d) != 0 )
    
    # Caso base: es el ideal total
    if esO(ideal,d):
        return {}

    # Caso recursivo: calcula el primer factor primo de la norma del
    # ideal, comprueba sus divisores y factoriza recursivamente con el
    # primero que lo divida.
    p = next(iter(factorint(normaIdeal(ideal,d))))
    for divisor in idealesDivisores(p,d):
        if divideIdeal(ideal,divisor,d):
            ideal = cocienteIdealPrimo(ideal,divisor,d)
            return dict(
                Counter({tuple(divisor):1}) +
                Counter(factoriza_id(ideal,d))
            )


def esprimo(ideal,d):
    """Comprueba si un ideal dado es primo en un dominio O(d).
    
    Dividimos por casos
    
        * si I es el ideal total, sabemos que no es primo.
        * si norma(I) es primo de Z, es primo.
        * si norma(I) es cuadrado de primo. El ideal es primo si y sólo si
          el p dado no ramifica. Un primo ramifica si y sólo si el polinomio
          irreducible de e, Irr(e), es reducible; lo que se comprueba observando
          si tiene raíces.
        * en otro caso, no es primo.
    """
    if isprime(normaIdeal(ideal,d)):
        # Los ideales de norma prima son primos
        return True
    else:
        p = sqrt(normaIdeal(ideal,d))
        if isprime(p):
            irr = irr_e(d)
            solucionesirr = [x for x in xrange(p) if irr(x)%p==0]
            ramifica = (len(solucionesirr) != 0)
            return not ramifica
        return False


def cocienteIdeal(I,J,d):
    """No sabemos calcular inversos de ideales en general, sólo sabemos
    calcular inversos de ideales principales generados por un primo;
    así que factorizaremos el ideal J para dividir por cada uno de sus
    factores a I.
    """
    resultado = I
    for k,v in factoriza_id(J,d).items():
        for x in xrange(v):
            resultado = cocienteIdealPrimo(resultado,k,d)
            
    return resultado


def cocienteIdealPrimo(I,u,d):
    """Calcula el ideal de un cociente por un ideal primo.
    """
    # Calcula el cociente dependiendo del tipo de divisor que
    # tengamos. Puede ser:
    #
    #  * un ideal principal generado por un primo. Con norma el
    #    cuadrado de un primo.
    #
    #  * un ideal generado por un primo y un elemento e-a. Con norma
    #    prima.
    #
    #assert( esprimo(u,d) )
    #assert( divideIdeal(I,u,d) )
    
    if not isprime(normaIdeal(u,d)):
        # El primo que determina el ideal es la raíz de la norma
        p = int(sqrt(normaIdeal(u,d)))
        return map(lambda i: i*Rational(1,p), I)
    else:
        p = normaIdeal(u,d)
        # Calcula la otra raíz del polinomio para buscar el inverso
        # por el que debe multiplicar para calcular el cociente.
        irr = irr_e(d)
        solucionesirr = [x for x in xrange(p) if irr(x)%p==0]
        for s in solucionesirr:
            a = [p, e(d)-s]
            if esO(map(lambda i: i*Rational(1,p),productodos(a,u,d)), d):
                return map(lambda i: i*Rational(1,p), productodos(I,a,d))










##
# PRÁCTICA 5. Número de clase
##

# Cálculo de la cota de Minkowski
def discriminante(d):
    """Calcula el discriminante del cuerpo."""
    if d%4 != 1:
        return 4*d
    else:
        return d


def minkowski(d):
    """Calcula la cota de Minkowski."""
    # Definición de s y t
    if d > 0:
        s = 2
        t = 0
    else:
        s = 0
        t = 1
    n = s + 2*t
    
    cota = (4/float(pi))**t * factorial(n)/(n**n) * sqrt(discriminante(d))
    return float(cota)


def ramifica(p,d):
    """Calcula si el primo ramifica en el cuerpo dado."""
    irr = irr_e(d)
    solucionesirr = [x for x in xrange(p) if irr(x)%p==0]
    return (len(solucionesirr) != 0)


def primosdebajode(x):
    """Calcula los primos por debajo de una cota."""
    return [x for x in sieve.primerange(0,x)]


def primosqueramificanbajolacota(d):
    """Devuelve los primos que ramifican debajo de la cota
    de Minkowski del cuerpo."""
    primos = primosdebajode(minkowski(d))
    return filter(lambda p: ramifica(p,d), primos)


def idealesGeneradores(d):
    """Los generadores son los ideales con la norma de los primos que
    ramifican.  Podemos coger sólo el primero porque si hubiera otro,
    serían inversos.
    """
    return map(lambda p: idealesDivisores(p,d)[0],
               primosqueramificanbajolacota(d))


def esPrinc(ideal,d):
    """Determina si el ideal dado es principal."""
    n = normaIdeal(ideal,d)
    l = chain(connorma(n,d), connorma(-n,d))
    return any(pertenece(u,ideal,d) for u in l)

def quitaTriviales(lista,d):
    """Quita los ideales triviales de una lista de generadores;
    es decir, aquellos que son principales."""
    return filter(lambda i: not esPrinc(i,d), lista)


def quitaInversos(lista,d):
    """Quita los inversos de una lista de generadores."""
    # Simplemente comprueba la multiplicación de todas las
    # parejas y va eliminando las inversas.
    i = 0
    while i < len(lista):
        j = i+1
        while j < len(lista):
            if esPrinc(productodos(lista[i],lista[j],d),d):
                print "El ideal",i,"ésimo es inverso del ideal",j,"ésimo."
                del lista[j]
            else:
                j = j+1
        i = i+1
    return lista

def quitaRelacionesTriviales(lista,d,n=2,debug=False):
    """Quita los generadores innecesarios de la lista de
    generadores. Un generador innecesario es aquel que puede
    escribirse como combinación de los demás. Esto es, aquel
    que cumple una relación teniendo coeficiente 1 en ella."""

    j = len(lista)-1
    while j >= 0:
        t = len(lista)

        # Prueba las relaciones posibles teniendo un 1 en la posición
        # dada por el índice j
        for relacionposible in sumasposibles(range(n+1), n, t-1):
            relacionposible.insert(j,1)

            if esPrinc(calculaRelacion(relacionposible,lista,d),d):
                if debug:
                    print "Relación encontrada",relacionposible,"retirado el",j,"ésimo generador"
                    
                del lista[j]
                break
        
        j = j-1
    return lista

def refinaGeneradores(lista,d,nmax=4):
    """Refina la lista de ideales generadores, reduciéndola en lo posible
    eliminando aquellos que cumplen relaciones triviales."""
    lista = quitaTriviales(lista,d)
    #lista = quitaInversos(lista,d)
    for n in xrange(nmax):
        quitaRelacionesTriviales(lista,d,n)
    return lista


def ordenIdeal(ideal, d):
    """Calcula el orden de un ideal en el grupo de clase del cuerpo."""
    n = 1
    i = ideal
    while not esPrinc(i,d):
        n = n+1
        i = productodos(i,ideal,d)
    return n


def ordenLista(lista, d):
    """Calcula el orden de cada elemento de una lista de generadores."""
    return map(lambda i: (i,ordenIdeal(i,d)), lista)


def matrizRelacionesOrden(ordenes):
    """Devuelve la matriz que genera el orden de los elementos."""
    relationlist = map(lambda (a,b): b, ordenes)
    relationmatrix = [
        [x if i==j else 0 for j in range(len(relationlist))] 
        for (i,x) in list(enumerate(relationlist))]

    return relationmatrix



def calculaRelacion(rel, generadores, d):
    """Calcula el resultado de multiplicar una relación con los
    generadores dados."""
    resultado = [1]
    
    for i in xrange(len(rel)):
        for _ in xrange(rel[i]):
            resultado = productodos(resultado, generadores[i], d)
    return resultado


def esNuevaRelacion(rel, matriz, generadores, d):
    """Comprueba si la relación dada es efectivamente una
    relación y si es además nueva. Esto es, no puede deducirse
    de las anteriores."""
    for relacion in matriz:
        if all([a >= b for (a,b) in zip(rel,relacion)]):
            return False
    
    if esPrinc(calculaRelacion(rel, generadores, d),d):
        return True
    
    return False

def sumasposibles(r,n,t):
    """Calcula las formas posibles de sumar n con t términos de la
    lista r."""
    return filter(lambda l: sum(l)==n, map(list,product(r,repeat=t)))

def sumasposibles2(r, n, t):
    """Calcula las formas posibles de sumar n con t términos de la
    lista r."""
    if t == 0:
        if n == 0:
            return [[]]
        return []
    for x in r:
        if x > n:
            break
        [suma + [x] for suma in sumasposibles2(r, n-x, t-1)]



def encuentraRelaciones(matriz, generadores, d, debug=False, n=3):
    """Busca las relaciones que pueden encontrarse con estos generadores,
    conociendo las que ya existen en la matriz."""
    t = len(generadores)
    alguncaso = True
    
    while alguncaso:
        if debug:
            print "Buscando relaciones con n =",n
        
        alguncaso = False
        for suma in sumasposibles(range(n+1),n,t):
            if esNuevaRelacion(suma, matriz, generadores, d):
                # Si encuentra una nueva relación, la añade a la matriz de
                # relaciones conocidas, simplifica y sigue buscando.
                matriz = [suma] + matriz
                print suma
                alguncaso = True
        n = n+1
    
    return matriz


def retiraFilasDeCeros(matriz):
    j=0
    while j < len(matriz):
        if all([x == 0 for x in matriz[j]]):
            del matriz[j]
        else:
            j = j+1
    return matriz

def reduceRelaciones(matriz):
    """Reduce una matriz de relaciones por filas."""
    n = len(matriz)
    if n == 0: return matriz
    m = len(matriz[0])
    
    for step in range(m):

        # mientras que la matriz no esté reducida
        while any([fila[step] != 0 for fila in matriz[(step+1):]]):

            # Ordena la matriz por valor absoluto, dejando ceros al final
            matriz[step:] = sorted(matriz[step:],
                key=lambda x: abs(x[step]) if x[step] != 0 else float("inf"))

            # Y simplifica por filas
            pivote = matriz[step][step]
            for i in range(step+1,n):
                coc = matriz[i][step] / pivote        
                for j in range(m):
                    matriz[i][j] = matriz[i][j] - matriz[step][j] * coc

    # Acaba retirando las innecesarias filas de ceros
    return retiraFilasDeCeros(matriz)



def numeroClase(d, debug=True):
    # Comprueba que el número sea libre de cuadrados.
    assert( libreDeCuadrados(d) )

    # Calcula la cota de minkowski y los ideales generadores
    cotaminkowski = minkowski(d)
    listaprimos = primosqueramificanbajolacota(d)
    ideales = idealesGeneradores(d)

    if debug:
        print "La cota de Minkowski para",d,"es",cotaminkowski
        print "Los primos que ramifican por debajo de esa cota son:"
        print listaprimos,"\n\n"
        print "Tenemos los ideales generadores"
        print ideales
    
    # Refina la lista de generadores
    lista = refinaGeneradores(ideales,d)
    
    if debug:
        print "Refinamos la lista de generadores retirando los inversos y nos quedamos con",lista,"\n\n"

    # Crea la matriz con el orden
    ordenes = ordenLista(lista,d)
    matriz = matrizRelacionesOrden(ordenes)

    if debug:
        print "El orden de los elementos nos da la siguiente matriz de relaciones"
        print np.matrix(matriz)

    # Busca más relaciones
    nuevamatriz = encuentraRelaciones(matriz,lista,d,debug=debug)

    # Reduce por filas la matriz hasta obtener los factores fundamentales
    final = reduceRelaciones(nuevamatriz)

    if debug:
        print "La matriz reducida es"
        print np.matrix(final)

    # Calcula el número de clase
    clase = 1
    if len(final)>0:
        for i in xrange(len(final[0])):
            clase = clase * abs(final[i][i])

    if debug:
        print "El número de clase es ",clase
        
    return clase
