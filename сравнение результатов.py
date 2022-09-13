from numpy import array, linspace, zeros, sin, trapz, cos, exp, tan
from numpy.linalg import solve
from sympy import diff, symbols, lambdify
from matplotlib.pylab import figure, plot, show, legend

def mkr(n,a,b,ya,yb):
    def f(x):
        return 3 * x - 4
    h = (b-a)/ n #шаг разностной сетки
    x = [i*h for i in range(1,n)] 
    A = zeros((n-1,n-1))
    B = zeros(n-1)
    A[0,0] = 8 * h**2 - 2; A[0,1] = 1 - 2 * h
    A[n-2,n-3] = 1 + 2 * h; A[n-2,n-2] = 8 * h**2 - 2
    B[0] = h**2 * f(x[0]) - ya
    B[n-2] = h**2 * f(x[n-2]) - yb
    for i in range(1,n-2):
        A[i,i-1] = 1 + 2 * h
        A[i,i] = 8 * h**2 - 2
        A[i,i+1] = 1 - 2 * h
        B[i] = h**2 * f(x[i])  
    y = solve(A,B) 
    X = zeros((n+1));X[0] = a;X[n] = b
    Y = zeros((n+1));Y[0] = ya;Y[n] = yb
    for i in range(1,n):
        X[i] = x[i-1]
        Y[i] = y[i-1]
    return X,Y    

def coll(n,a,b,ya,yb,nx):
    h = (b - a) / n #шаг разностной сетки
    xi = [i*h for i in range(1,n)] 
    A = zeros((n-1,n-1))
    B = zeros(n-1)
    for i in range(n-1):
        xx = xi[i]
        for j in range(1, n):
            l, m = symbols('l, m')
            phi = l ** m * (1 - l); dphi = diff(phi, l); ddphi = diff(dphi, l)
            phi = lambdify((l, m), phi)
            dphi = lambdify((l, m), dphi)
            ddphi = lambdify((l, m), ddphi)
            A[i,j-1] = ddphi(xx, j) - 4 * dphi(xx, j) + 8 * phi(xx, j)
        B[i] =  -xx - 6
    ab = solve(A,B)   
    
    x = linspace(a,b,nx)
    y = zeros(nx);y[0] = ya; y[nx-1] = yb
    for i in range(1,nx-1):
        xx = x[i] #;print(xx)
        y[i] = 0.5*(xx + 1)
        for j in range(1,n):
           y[i] += ab[j-1] * xx ** j * (1 - xx)
    return x,y

def bubgal(n,a,b,ya,yb,nx):
    A = zeros((n-1,n-1))
    B = zeros((n-1))
    xi = linspace(a,b,n+1)
    for i in range(1, n):
      x = xi[i]
      p =  xi ** i * (1 - xi)
      l, m = symbols('l, m')
      phi = l ** m * (1 - l)
      dphi = diff(phi, l)
      ddphi = diff(dphi, l)
      phi = lambdify((l, m), phi)
      dphi = lambdify((l, m), dphi)
      ddphi = lambdify((l, m), ddphi)
      for j in range(1, n):
        A[i-1,j-1] = trapz(p * (ddphi(x,j)-4*dphi(x,j)+8*phi(x,j)), xi)
      B[i-1] = trapz(p * (-x - 6), xi)
    ab = solve(A,B)
    x = linspace(a,b,nx)
    y = zeros((nx));y[0] = ya; y[nx-1] = yb
    for i in range(1,nx-1):
        xx = x[i]
        #print(xx)
        y[i] = 0.5*(xx + 1)
        for j in range(1,n):
            y[i] += ab[j-1] * xx ** j * (1 - xx)
    return x,y

def mke(n,a,b,ya,yb):
    def integ(a,b,m):
        x = linspace(a,b)
        if m: l = x - a
        else: l = b - x
        return trapz(l/(b-a)*(3 * x - 4), x)
    h = (b - a) / n  
    E = n + 1 #количество узлов
    A = zeros((E,E))
    B = zeros((E))
    x = [i*h for i in range(n+1)]
    Ak =  1 / h *array([[1, -1], [-1, 1]]) + 2 * array([[-1, 1], [-1, 1]])\
            - 8 * h * array([[1/3, 1/6], [1/6, 1/3]])
    for i in range(E-1):
        A[i:i+2,i:i+2] += Ak
    for i in range(E-1):
        x_k = x[i]; x_k1 = x[i+1]
        Bk = [0, 0]
        Bk[0] = integ(x_k, x_k1, 0)
        Bk[1] = integ(x_k, x_k1, 1)
        B[i:i+2] += Bk
    Aa = A[1:n,1:n] 
    Bb = -B[1:n].copy()  
    Bb[0] -= A[1,0]*ya; Bb[n-2] -= A[n-1,n]*yb
    ab = solve(Aa,Bb)
    X = zeros((n+1));X[0] = a;X[n] = b
    Y = zeros((n+1));Y[0] = ya;Y[n] = yb
    for i in range(1,n):
        X[i] = x[i]
        Y[i] = ab[i-1]
    return X,Y

n = 10
nx = 1000
a, b = 0, 1
ya, yb = 0.5, 1
rgb = ('r','g','b', 'y')
figure()
x0, y0 = mkr(n,a,b,ya,yb)
plot(x0,y0, color = rgb[0], label = f'МКР N = {n}')
x1, y1 = coll(n,a,b,ya,yb,nx)
plot(x1,y1, color = rgb[1], label = f'Метод коллокаций N = {n}')
x2, y2 = bubgal(n,a,b,ya,yb,nx)
plot(x2,y2, '--', color = rgb[2], label = f'Метод Бубнова-Галеркина N = {n}')
x3, y3 = mke(n,a,b,ya,yb)
plot(x3,y3, color = rgb[3], label = f'МКЭ N = {n}')

xt = linspace(a,b,nx+1)
c1 = 13/16
c2 = 15/16 * exp(-2) * 1 / sin(2) - 13/16 * 1 / tan(2)
yt = exp(2*xt)*(c1*cos(2*xt)+c2*sin(2*xt)) + 3/8*xt - 5/16
plot(xt,yt,'--', linewidth = 2, color = 'black', label = 'Аналитическое решение')     
legend()
show() 