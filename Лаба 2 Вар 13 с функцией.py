from numpy import cos, trapz, linspace, zeros,size
from matplotlib import pylab as plt

def Scalar(x,f,g):                 #функция скалярного умножения
    return trapz(f*g,x)

def Norm(x,f):                     #функция подсчета нормы
    return trapz(f**2,x)**0.5

def grsh(PHI,xs,N):                #функция Грама - Шмидта
    PSI = zeros((size(xs),N))
    PSI[:,0] = PHI[:,0]/Norm(xs,PHI[:,0])
    for i in range(1, N):
        FScalar = 0
        for j in range(0, i):
            FScalar -= Scalar(xs,PHI[:, i], PSI[:, j])*PSI[:, j]
        F = PHI[:,i] + FScalar
        PSI[:,i] = F / Norm(xs, F)
    return PSI
            
#Ввод функций
x = linspace(0,1)                   #задание x
listPhi = [1,cos(x),x,x**2,x**3-1]  #Создание списка с функциями phi
N = len(listPhi)
PHI = zeros((size(x),N))
for k in range(N):                  #добавление в массив функций из списка
    PHI[:,k] = listPhi[k]
            
PSI = grsh(PHI,x,N)                 #подсчет PSI через функцию Грама-Шмидта

#вывод графиков
plt.figure()
for i in range(N):
    plt.plot(x,PHI[:, i])
plt.figure()
for i in range(N):
    plt.plot(x,PSI[:, i])
    
print(Scalar(x,PSI[:,2],PSI[:,3]))