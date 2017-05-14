#! /bin/python2
# -*- encoding: utf-8 -*-

from sympy import *
from collections import Counter


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
    cuerpo cuadrático O(d). Devuelve un natural en la salida.

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
    return ask(Q.integer(norma(a, d))) and ask(Q.integer(traza(a, d)))


def xy(a,d):
    """Escribe las coordenadas del número del cuerpo algebraico en parte
    racional y coeficiente del radical.

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
    solucionespotenciales = ((sqrt(n+d*y**2), y) for y in xrange(cotainf,cotasup+1))
    soluciones = [
        [(x,y),(x,-y),(-x,y),(-x,-y)]
        for (x,y) in solucionespotenciales if ask(Q.integer(x))
    ]
    
    return list(set([s for sol in soluciones for s in sol]))


def eqpell_neg(n,d):
    """Resuelve la ecuación de Pell para d<0.

    """
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
    solucionespotenciales = ((sqrt(n+d*y**2), y) for y in xrange(cota))
    soluciones = [[(x,y),(x,-y),(-x,y),(-x,-y)] for (x,y) in solucionespotenciales if ask(Q.integer(x))]

    # Devolvemos las soluciones sin repetición
    return list(set([s for sol in soluciones for s in sol]))


def eqpell(n,d):
    """Resuelve la ecuación de Pell. En el caso de que d sea positivo sólo
    devolverá las soluciones generadoras.

    """
    if n == 1 and d > 0:
        return pell(d)
    
    if d <= 0:
        return eqpell_neg(n,d)
    else:
        return generalpell(n,d)


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
    
    # Comprueba si la raíz de un número es primo. Si lo es, intenta
    # comprobar si ramifica buscando soluciones.
    raiz = sqrt(n)
    if ask(Q.integer(raiz)) and isprime(raiz):
        return len(connorma(raiz,d) + connorma(-raiz,d)) == 0
    
    return False


def factoriza(a,d):
    """Factoriza el elemento en el cuerpo O(d).

    """
    assert( type(a) != list )
    
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
    connormadada = connorma(factor, d) + connorma(-factor,d)
    if len(connormadada) != 0:
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
    if d % 4 == 1:
        return Rational(1,2)+Rational(1,2)*sqrt(d) 
    else:
        return sqrt(d)

def irr_e(d):
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
    reducida = LR(matrizrel(ideal,d))
    
    # Caso particular de longitud 1 o 0.
    if len(reducida) < 2:
        return reducida
    
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
    return simplificaIdeal([i*j for i in I for j in J],d)


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
    assert( esprimo(u,d) )
    assert( divideIdeal(I,u,d) )
    
    if not isprime(normaIdeal(u,d)):
        p = u[0]
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
