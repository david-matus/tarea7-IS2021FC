import numpy as np
import os
import matplotlib.pyplot as plt
import time

def CrearDirectorio(folder):
    '''
    Crea una carpeta con el nombre dado en la ubicación donde se está ejecutando el programa.
    Si existe una carpeta con el mismo nombre, crea una con la forma nombre (1), etc.
    Entradas:
    folder: str: nombre deseado de la carpeta
    Salidas:
    subfolder: str: nombre con el que fue creada la carpeta
    '''
    try: #primero se intenta guardar con un nombre sencillo
        os.mkdir(folder)
        subfolder = folder
    except:
        filepend = True
        trycount = 1
        while filepend:
            try: #se añade una cola al nombre de carpeta si la carpeta ya existe
                os.mkdir(folder+' ('+str(trycount)+')')
                filepend = False
                subfolder = folder+' ('+str(trycount)+')'
            except:
                trycount += 1
    return subfolder #el nombre de la carpeta creada

def GraficarRuta(camino,ciudadesXY,generación,subdirectorio):
    '''
    Crea una gráfica con matplotlib de la trayectoria seguida por la secuencia de un cromosoma.
    Entradas:
    camino: 1xn matrix: camino a graficar
    ciudadesXY: 2xn matrix: matriz [X,Y] que contiene las coordenadas de cada ciudad
    generación: int: ciclo principal de la simulación
    subdirectorio: str: nombre de la carpeta donde se guardará la gráfica
    Salidas:
    Esta función no retorna ningún valor, pero crea un archivo png en el directorio indicado.
    '''
    secuenciasub = np.copy(camino)
    secuencia = np.append(secuenciasub,secuenciasub[0])
    totalGenes = len(secuencia)
    for gen in range(totalGenes-1):
        ciudad1 = secuencia[gen]
        ciudad2 = secuencia[gen+1]
        x1, x2 = ciudadesXY[0][ciudad1], ciudadesXY[0][ciudad2]
        y1, y2 = ciudadesXY[1][ciudad1], ciudadesXY[1][ciudad2]
        xstep = abs(x2-x1)/1000
        ystep = abs(y2-y1)/1000
        if x1<x2:
            x = np.arange(x1,x2,xstep)
            if y1<y2:
                y = np.arange(y1,y2,ystep)
            elif y1>y2:
                y = np.flip(np.arange(y2,y1,ystep))
            else:
                y = np.ones(1000)
                y[:]*=y1
        elif x1>x2:
            x = np.arange(x2,x1,xstep)
            if y1<y2:
                y = np.flip(np.arange(y1,y2,ystep))
            elif y1>y2:
                y = np.arange(y2,y1,ystep)
            else:
                y = np.ones(1000)
                y[:] *= y1
        else:
            x = np.ones(1000)
            x[:] *= x1
            if y1<y2:
                y = np.arange(y1,y2,ystep)
            else:
                y = np.flip(np.arange(y2,y1,ystep))
        plt.plot(x[0:1000],y[0:1000],'b')
    plt.plot(ciudadesXY[0],ciudadesXY[1],'r*')
    genes = len(camino)
    plt.plot(ciudadesXY[0][(camino[0])],ciudadesXY[1][(camino[0])],'g*',label='inicio')
    plt.plot(ciudadesXY[0][(camino[genes-1])],ciudadesXY[1][(camino[genes-1])],'y*',label='fin')
    plt.title('Mejor ruta en gen: '+str(generación))
    endPath = subdirectorio+'/'+'Ruta en gen '+str(generación)+'.png'
    plt.savefig(endPath)
    plt.close('all')
    return

def LeerArchivo(nombreTxt):
    '''
    Obtiene los puntos en XY de las ciudades a utilizar a partir de una matriz almacenada en un archivo
    txt. Los retorna como una matriz [X,Y]
    Entradas:
    nombreTxt: str: nombre del archivo .txt que se va a leer
    Salidas:
    [X,Y]: 2xn matrix: matriz [X,Y] que contiene las coordenadas de cada ciudad
    '''
    datosTxt = open(nombreTxt,'rt')
    listaDeStr = (datosTxt.read()).split('\n')
    for m in range(len(listaDeStr)):
        string = listaDeStr[m]
        listaDeStr[m] = string.replace(']','').replace('[','').replace(' ','')
    for l in range(len(listaDeStr)):
        string = listaDeStr[l]
        listaDeStr[l] = string[0:(len(string)-1)]
    listaDeStr = listaDeStr[0:(len(listaDeStr)-1)]
    X = []
    Y = []
    for x in range(len(listaDeStr)):
        X.append(float(listaDeStr[x].split(',')[0]))
        Y.append(float(listaDeStr[x].split(',')[1]))
    return [X,Y]

def EscribirArchivo(matriz,nombreTxt):
    '''
    Escribe un archivo txt de la matriz que se indica en el mismo directorio que este script.
    Entradas:
    matriz: numpy array: array que se debe almacenar en el archivo txt
    nombreTxt: str: nombre del archivo .txt que se va a almacenar
    Salidas:
    Esta función no retorna ningún valor, pero crea un archivo txt en el directorio del script.
    '''
    with open(nombreTxt,'w') as archivo:
        archivo.write(str(matriz))
    return

def DistanciaEntreCiudades(nodo1, nodo2, coordCiudades):
    '''
    Calcula la distancia euclideana entre dos nodos del plano XY.
    Entradas:
    nodo1: int: índice del primer nodo
    nodo2: int: índice del segundo nodo
    coordCiudades: 2xn matrix: matriz [X,Y] con las coordenadas de los nodos
    Salidas:
    longitudTrayectoria: float: distancia euclideana entre los nodos
    '''
    dx = coordCiudades[0][nodo1]-coordCiudades[0][nodo2]
    dy = coordCiudades[1][nodo1]-coordCiudades[1][nodo2]
    longitudTrayectoria = np.sqrt(dx**2+dy**2)
    return longitudTrayectoria

def ObtenerVisibilidad(coordCiudades):
    '''
    Calcula la matriz de visibilidad entre los diferentes nodos.
    Entradas:
    coordCiudades: 2xn matrix: matriz [X,Y] con las coordenadas de los nodos
    Salidas:
    arregloVisibilidad: nxn matrix: matriz de visibilidad
    '''
    nciudades = len(coordCiudades[0])
    arregloVisibilidad = np.zeros((nciudades,nciudades))
    for iciudad1 in range(nciudades):
        for iciudad2 in range(nciudades):
            if iciudad1 != iciudad2:
                d = DistanciaEntreCiudades(iciudad1, iciudad2, coordCiudades)
                #visibilidad es inverso de distancia euclideana
                arregloVisibilidad[iciudad1][iciudad2] = 1/d
    return arregloVisibilidad

def ObtenerProbabilidad(nivelFeromonas, visibilidad, alpha, beta, listaTabú, iNodo, jNodo):
    '''
    Calcula la probabilidad de pasar de nodo i a nodo j.
    Entradas:
    nivelFeromonas: nxn matrix: matriz con las feromonas entre nodos
    visibilidad: nxn matrix: matriz con la visibilidad entre nodos
    alpha, beta: float, float: parámetros de la simulación
    listaTabú: lista de nodos visitados
    iNodo: índice del nodo i
    jNodo: índice del nodo j
    Salidas:
    probabilidadNodo: float: probabilidad de pasar de nodo i a nodo j
    '''
    nnodos = len(nivelFeromonas[0])
    denominador = 0
    for t in range(nnodos):
        if t not in listaTabú:
            denominador += (nivelFeromonas[iNodo][t]**alpha)*(visibilidad[iNodo][t]**beta)
    probabilidadNodo = (nivelFeromonas[iNodo][jNodo]**alpha)*(visibilidad[iNodo][jNodo]**beta)/denominador
    return probabilidadNodo

def GenerarCamino(nivelFeromonas, visibilidad, alpha, beta):
    '''
    Calcula la probabilidad de pasar de nodo i a nodo j.
    Entradas:
    nivelFeromonas: nxn matrix: matriz con las feromonas entre nodos
    visibilidad: nxn matrix: matriz con la visibilidad entre nodos
    alpha, beta: float, float: parámetros de la simulación
    Salidas:
    caminoGenerado: 1xn list: lista de índices de la ruta
    '''
    nnodos = len(nivelFeromonas[0])
    nodo0 = np.random.randint(0,nnodos)
    listaTabú = [nodo0]
    caminoGenerado = [nodo0]
    for x in range(nnodos-1): #cada ciclo es un nuevo punto en el camino
        sumaParcial = 0
        choose = np.random.random()
        for y in range(nnodos): #cada ciclo es una nueva evaluación de P
            if y not in listaTabú:
                sumaParcial += ObtenerProbabilidad(nivelFeromonas, visibilidad, alpha, beta, listaTabú, caminoGenerado[x], y)
                if choose>(1-sumaParcial):
                    caminoGenerado.append(y)
                    listaTabú.append(y)
                    break
    return caminoGenerado

def CálculoDeltaTau(colecciónCaminos, ciudadesXY):
    '''
    Calcula el cambio de la matriz de feromonas debido al depóstito de feromonas por las hormigas.
    Entradas:
    colecciónCaminos: mxn matrix: matriz con las rutas de todas las hormigas
    ciudadesXY: 2xn matrix: matriz [X,Y] con las coordenadas de las ciudades
    Salidas:
    deltaTau: nxn matrix: matriz con el cambio de tau ij
    '''
    nnodos = len(ciudadesXY[0])
    nhormigas = len(colecciónCaminos)
    deltaTau = np.zeros((nnodos,nnodos))
    for k in range(nhormigas): #cada hormiga deposita feromonas
        for m in range(nnodos): #las feromonas se depositan en cada fragmento ij de la ruta
            i = colecciónCaminos[k][m-1]
            j = colecciónCaminos[k][m]
            distancia = DistanciaEntreCiudades(i, j, ciudadesXY)
            deltaTau[i][j] += 1/distancia
            deltaTau[j][i] += 1/distancia #la matriz de feromonas debe ser simétrica                
    return deltaTau

def ActualizarNivelFeromonas(nivelFeromonas, deltaNivelFeromonas, rho):
    '''
    Establece los nuevos valores de la matriz de feromonas ij.
    Entradas:
    nivelFeromonas: nxn matrix: matriz con las feromonas entre nodos
    deltaNivelFeromonas: nxn matrix: matriz con los valores de delta tau para cada ij
    rho: float: parámetro de la simulación
    Salidas:
    nivelFeromonas: nxn matrix: nueva matriz con las feromonas entre nodos
    '''
    for x in range(len(nivelFeromonas)):
        for y in range(len(nivelFeromonas)):
            nivelFeromonas[x][y] = (1-rho)*nivelFeromonas[x][y]+deltaNivelFeromonas[x][y]
    return nivelFeromonas

def EvaluarPoblación(población,ciudadesXY):
    '''
    Calcula el valor de evaluación de cada uno de las rutas en la población. La evaluación
    es el inverso de la distancia total recorrida por la secuencia de la ruta.
    Entradas:
    población: nxm matrix: matriz que contiene las rutas
    ciudadesXY: 2xn matrix: matriz [X,Y] que contiene las coordenadas de cada ciudad
    Salidas:
    puntuaciones: 1xn numpy array: lista de puntuaciones de cada ruta
    '''
    Totalcromosomas = len(población)
    Totalgenes = len(población[0])
    puntuaciones = np.zeros(Totalcromosomas)
    for icromosoma in range(Totalcromosomas):
        puntuacionActual = 0
        cromosomaMod = np.copy(población[icromosoma])
        cromosomaActual = np.append(cromosomaMod,cromosomaMod[0])
        for igen in range(Totalgenes):
            a = cromosomaActual[igen]
            b = cromosomaActual[igen+1]
            x1, y1 = ciudadesXY[0][a], ciudadesXY[1][a]
            x2, y2 = ciudadesXY[0][b], ciudadesXY[1][b]
            puntuacionActual += np.sqrt((x1-x2)**2+(y1-y2)**2)
        puntuaciones[icromosoma] = 1/puntuacionActual
    return puntuaciones

def EncontrarMaxPop(población,puntuaciones):
    '''
    Calcula el índice dentro de la población de la ruta con mayor puntuación.
    Entradas:
    población: nxm matrix: matriz que contiene las rutas
    puntuaciones: 1xn numpy array: lista de puntuaciones de cada ruta
    Salidas:
    imax: int: índice dentro de la población de la mejor ruta
    '''
    imax = 0
    for i in range(len(población)): #se encuentra el mejor cromosoma
        if puntuaciones[i]>puntuaciones[imax]:
            imax = i
    return imax

def Optimización():
    np.random.seed(12321)
    #Parámetros
    nHormigas = 20
    nIteraciones = 1000
    alpha = 1.0
    beta = 5.0
    rho = 0.5
    Lobj = 120
    nombreTxt = 'CoordenadasCiudades.txt'
    #Inicializar variables
    ciudadesXY = LeerArchivo(nombreTxt)
    nnodos = len(ciudadesXY[0])
    tau = np.ones((nnodos,nnodos)) #nivel de feromonas
    eta = ObtenerVisibilidad(ciudadesXY) #visibilidad
    subdirectorio = CrearDirectorio('Gráficas')
    subdirectorioR = CrearDirectorio('Resultados')
    puntuacionObjetivo = 1/Lobj
    mejorPuntuacion = 0
    t = 0
    lenMejorRuta = 0
    mejoresLongitudes = []
    LongitudesProm = []
    rutaOptima = []
    while t<nIteraciones and mejorPuntuacion<puntuacionObjetivo:
    #for t in range(nIteraciones): #cada ciclo es un recorrido por todas las hormigas
        colecciónRutas = []
        for k in range(nHormigas): #cada ciclo es el recorrido de una hormiga
            #mover a la hormiga según P
            rutaK = GenerarCamino(tau,eta,alpha,beta)
            colecciónRutas.append(rutaK)
        #actualizar tau usando dtau
        dtau = CálculoDeltaTau(colecciónRutas, ciudadesXY)
        tau = ActualizarNivelFeromonas(tau, dtau, rho)
        #Averiguar si alguna ruta es la mejor hasta el momento
        puntuaciones = EvaluarPoblación(colecciónRutas,ciudadesXY)
        iMejorActual = EncontrarMaxPop(colecciónRutas,puntuaciones)
        lenMejorActual = 1/puntuaciones[iMejorActual]
        if mejorPuntuacion < puntuaciones[iMejorActual]:
            #Graficar mejor ruta
            GraficarRuta(colecciónRutas[iMejorActual],ciudadesXY,t,subdirectorio)
            lenMejorRuta = lenMejorActual
            mejorPuntuacion = puntuaciones[iMejorActual]
            rutaOptima = colecciónRutas[iMejorActual]
            localtime = time.asctime( time.localtime(time.time()) )
            print('Tiempo de nueva mejor ruta:',localtime)
        mejoresLongitudes.append(lenMejorRuta)
        LongitudesProm = np.append(LongitudesProm,1/np.average(puntuaciones))
        t+=1
    generaciónvector = range(t)
    EscribirArchivo(np.array(rutaOptima),subdirectorioR+'/'+'caminoMásCorto_OSH.txt') #se guarda la ruta más corta
    GraficarRuta(rutaOptima,ciudadesXY,t,subdirectorio) #se grafica la ruta más corta
    print('Longitud mínima: '+str(mejoresLongitudes[len(mejoresLongitudes)-1]))
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(generaciónvector, mejoresLongitudes, 'b',label = 'Longitud más corta')
    ax.plot(generaciónvector, LongitudesProm, 'r',label = 'Longitud promedio')
    ax.set_title('Longitud vs generación')
    ax.legend()
    fig.savefig(subdirectorioR+'/'+'LvsGen OSH'+'.png')
    plt.close('all')
    return

localtime = time.asctime( time.localtime(time.time()) )
print('Tiempo de inicio :',localtime)
Optimización()
