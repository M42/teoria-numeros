#! /bin/python2
# -*- encoding: utf-8 -*-

from sympy import *
from FDEMRG import *


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
    generadores = [x for par in map(lambda a: [a, expand(e(d)* a)], ideal) for x in par]
    return map(lambda a: list(ab(a,d)), generadores)


def LR(matriz):
    """Calcula la forma reducida asociada a una matriz. La forma reducida
    será la matriz equivalente triangular inferior.

    """
    
    # Ordena la matriz por valor absoluto, dejando los ceros al final.
    matriz = sorted(matriz,
                    key=lambda x: abs(x[0]) if x[0] != 0 else
                    float("inf"))
    
    # Si ha simplificado ya la primera fila, simplifica sólo la
    # segunda fila.
    if (len(matriz) < 2 or matriz[1][0] == 0):
        if (len(matriz) < 3 or matriz[2][1] == 0):
            return matriz
        
        matriz = [matriz[0]] + sorted(matriz[1:],
                                      key=lambda x: abs(x[1]) if x[1] != 0 else
                                      float("inf"))
        a = matriz[1][1]
        for i in range(2,len(matriz)):
            matriz[i][1] = matriz[i][1] % a
        
        return LR(matriz)
        
    
    # Si no la ha simplificado, reduce desde el pivote.
    a = matriz[0][0]
    b = matriz[0][1]
    for i in range(1,len(matriz)):
        coc = matriz[i][0] / a
        matriz[i][0] = matriz[i][0] - a * coc
        matriz[i][1] = matriz[i][1] - b * coc
    
    return LR(matriz)


def simplificaIdeal(ideal, d):
    """Simplifica un ideal para escribirlo como generado por, a lo sumo,
    dos generadores.
    
    """
    reducida = LR(matrizrel(ideal,d))
    
    # Caso particular de longitud 1 o 0.
    if len(reducida) < 2:
        return reducida
    
    return [reducida[0][0] + reducida[0][1]*e(d), reducida[1][1]]


def normaIdeal(ideal, d):
    """Calcula la norma de un ideal. Viene dado como lista de
    generadores.

    """
    # Reduce la matriz de relatores del ideal. La norma se obtiene
    # como el producto de la primera subdiagonal.
    reducida = LR(matrizrel(ideal, d))
    return reducida[0][0] * reducida[1][1]
    
    
def esO(ideal, d):
    """Comprueba si un ideal es el total. Esto es, si es igual a O. El
    ideal viene dado como lista de generadores.
    
    """
    # Un ideal es el total si y sólo si tiene norma 1.
    return normaIdeal(ideal, d) == 1


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
        return false
    
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
    return simplificaIdeal([i*j for i in I for j in J],d)


def producto(listaideales):
    """Calcula el producto de una lista de ideales. Devuelve dos
    generadores del ideal producto.

    """
    return reduce(productodos, listaideales, [(1,0),(0,1)])


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
    return true


def cocienteIdeal(I,J,d):
    """No sabemos calcular inversos de ideales en general, sólo sabemos
    calcular inversos de ideales principales generados por un primo;
    así que factorizaremos el ideal J para dividir por cada uno de sus
    factores a I.

    """
    resultado = I
    for k,v in factoriza_id(J).items():
        for x in xrange(v):
            resultado = cocienteIdealPrimo(resultado,x,d)
            
    return resultado

def cocienteIdealPrimo(I,u,d):
    """Calcula el ideal de un cociente por un ideal primo.
    
    """
    # Calcula el cociente dependiendo del tipo de divisor que
    # tengamos. Puede ser:
    #  * un ideal principal generado por un primo. Con norma el
    #    cuadrado de un primo.
    #  * un ideal generado por un primo y un elemento e-a. Con norma
    #    prima.
    
    assert( esprimo(u,d) )
    if not isprime(normaIdeal(u,d)):
        p = u[0]
        return map(lambda i: i*Rational(1,p), I)
    else:
        p = normaIdeal(u,d)
        return map(lambda i: i*Rational(1,p), productodos(I,u,d))
