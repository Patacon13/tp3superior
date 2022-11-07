# Polinomio interpolación
# Diferencias Divididas de Newton
# Tarea: Verificar tamaño de vectores,
#        verificar puntos equidistantes en x
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

listaPolinomios = []

# PROCEDIMIENTO

# Tabla de Diferencias Divididas Avanzadas

def creaTabla(xi, fi):

    titulo = ['i   ','xi  ','fi  ']
    n = len(xi)
    ki = np.arange(0,n,1)
    tabla = np.concatenate(([ki],[xi],[fi]),axis=0)
    tabla = np.transpose(tabla)

    # diferencias divididas vacia
    dfinita = np.zeros(shape=(n,n),dtype=float)
    tabla = np.concatenate((tabla,dfinita), axis=1)

    return tabla, dfinita

# Calcula tabla, inicia en columna 3
def tablaDiferenciasDivididas(tabla, xi, fi):
    [n,m] = np.shape(tabla)
    diagonal = n-1
    j = 3
    while (j < m):
        # Añade título para cada columna
        #titulo.append('F['+str(j-2)+']')

        # cada fila de columna
        i = 0
        paso = j-2 # inicia en 1
        while (i < diagonal):
            denominador = (xi[i+paso]-xi[i])
            numerador = tabla[i+1,j-1]-tabla[i,j-1]
            if denominador == 0:
                xSymbol = sym.Symbol('x')
                polinomioDerivado = listaPolinomios[-1].diff(xSymbol)
                print("DERIVADA ", polinomioDerivado.evalf(subs={xSymbol:xi[0]}))

                tabla[i,j] = polinomioDerivado.evalf(subs={xSymbol:xi[0]})
            else:
                tabla[i,j] = numerador/denominador
            i = i+1
        diagonal = diagonal - 1
        j = j+1

        #print(tabla)

    return tabla

# POLINOMIO con diferencias Divididas
# caso: puntos equidistantes en eje x
def diferenciaDividida(tabla):
    dDividida = tabla[0,3:]
    return dDividida

def calculaPolinomio(dfinita, dDividida, xi, fi):
    n = len(dfinita)
    # expresión del polinomio con Sympy
    x = sym.Symbol('x')
    polinomio = fi[0]
    for j in range(1,n,1):
        factor = dDividida[j-1]
        termino = 1
        for k in range(0,j,1):
            termino = termino*(x-xi[k])
        polinomio = polinomio + termino*factor
    
    return polinomio.expand()

def polinomioTramo(xi, fi):
    tabla, dfinita = creaTabla(xi, fi)
    tablaDiferenciasDivididas(tabla, xi, fi)
    dDividida = diferenciaDividida(tabla)
    polinomio = calculaPolinomio(dfinita, dDividida, xi, fi)
    listaPolinomios.append(polinomio)

def cortaDeATramos(ptosX, ptosY, listaCortes):
    intervalosX = []
    intervalosY = []
    indice = 0
    for corte in listaCortes:
        esteIntervaloX = []
        esteIntervaloY = []
        print(indice)
        if indice > 0: 
            esteIntervaloX.append(ptosX[indice])
            esteIntervaloY.append(ptosY[indice])
            
        while corte >= ptosX[indice]:
            esteIntervaloX.append(ptosX[indice])
            esteIntervaloY.append(ptosY[indice])
            indice += 1
        intervalosX.append(esteIntervaloX)
        intervalosY.append(esteIntervaloY)

        indice -= 1
    esteIntervaloX = []
    esteIntervaloY = []

    indiceAux = 0

    while(indice<len(ptosX)):
        if indiceAux == 0:
            esteIntervaloX.append(ptosX[indice])
            esteIntervaloY.append(ptosY[indice])
            indiceAux += 1
        esteIntervaloX.append(ptosX[indice])
        esteIntervaloY.append(ptosY[indice])
        indice+=1
    
    intervalosX.append(esteIntervaloX)
    intervalosY.append(esteIntervaloY)
    
    return intervalosX, intervalosY

# INGRESO , Datos de prueba
xi = [0.00, 0.40, 0.80, 1.00, 1.15, 1.30, 1.50, 1.70, 1.90, 2.00, 2.10, 2.30, 2.40, 2.50, 2.60, 2.70, 3.00, 3.30, 3.60, 4.00, 4.50, 5.00, 5.50, 6.00]
fi = [-70.00, -70.00, -69.72, -65.78, -56.94, -48.28, -34.49, -15.21, 10.96, 29.44, 39.64, 14.19, -16.24, -45.10, -65.76, -78.98, -87.38, -84.70, -80.08, -75.12, -71.00, -70.00, -70.00, -70.00]

listaCortes = [0.8, 1.3, 1.9, 2.4, 2.6, 3.0, 4.5]

intervalosX, intervalosY = cortaDeATramos(xi, fi, listaCortes)

for i in range(0, len(intervalosX)):
    print(i)
    print(intervalosX[i])
    print(intervalosY[i])
    polinomioTramo(intervalosX[i], intervalosY[i])


#px = sym.lambdify(x,listaPolinomios[-1])

# Puntos para la gráfica
#muestras = 101
#a = np.min(xi)
#b = np.max(xi)
#pxi = np.linspace(a,b,muestras)
#pfi = px(pxi)

# SALIDA
np.set_printoptions(precision = 4)
print('Tabla Diferencia Dividida')
#print([titulo])
#print(tabla)
#print('dDividida: ')
#print(dDividida)
#print('polinomio: ')
#print(polinomio)
print('polinomio simplificado: ' )
print(listaPolinomios)

# Gráfica
plt.title("newton_interpolation")
plt.plot (xi, fi, 's') #El punto azul representa el valor original
xInicial = 0
xFinal = listaCortes[0]
for i in range(len(listaPolinomios)):
    '''
    x = np.linspace(0, 6)
    y = np.piecewise(x, [(x > 0) & (x <= 1.5), (x >= 1.5) & (x <= 2.5), (x >= 2.5) & (x <= 4.5), (x >= 4.5) & (x <= 6)], [])
    print(y)
    '''
    x = np.arange(xInicial, xFinal, 0.001)
    if (i + 1 < len(listaCortes)):
        xFinal = listaCortes[i + 1]
        xInicial = listaCortes[i]
    else:
        xFinal = intervalosX[-1][-1] 
        xInicial = listaCortes[-1]
    y = [listaPolinomios[i].evalf(subs={'x':val}) for val in x]
    plt.plot (x, y, 'r') # Curva de interpolación
    plt.xlabel('x')  
    plt.ylabel('y')  
    plt.legend (loc = 4) # Especifique la posición de la leyenda

plt.show()
'''
plt.plot(xi,fi,'o', label = 'Puntos')
##for i in range(0,n,1):
##    plt.axvline(xi[i],ls='--', color='yellow')
plt.plot(pxi,pfi, label = 'Polinomio')
plt.legend()
plt.xlabel('xi')
plt.ylabel('fi')
plt.title('Diferencias Divididas - Newton')
plt.show()
'''
