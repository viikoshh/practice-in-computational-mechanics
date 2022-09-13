from numpy import zeros, sin, cos, linspace, exp, tan
from numpy.linalg import solve
from matplotlib.pylab import figure, plot, show, legend

def f(x):
    return 3 * x - 4

def mkr(n,a,b,ya,yb):
    h = (b-a)/ n #шаг разностной сетки
    x = [i*h for i in range(1,n)] #;print(x)
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
    #print(A);print(B)   
    y = solve(A,B) #;print(y)
    X = zeros((n+1));X[0] = a;X[n] = b
    Y = zeros((n+1));Y[0] = ya;Y[n] = yb
    for i in range(1,n):
        X[i] = x[i-1]
        Y[i] = y[i-1]
    #print(X, Y)
    return X,Y

sn = (4, 6, 50, 1000) 
a, b = 0, 1
ya, yb = 0.5, 1
rgb = ('r','g','b', 'y')
#line = ('-','--',':')
#linew = (1, 2, 4)
figure()
for k in range(len(sn)):
    n = sn[k]
    x, y = mkr(n,a,b,ya,yb)
    plot(x,y, color = rgb[k], label = f'N = {n}')
#точное решение
nx = 10000
xt = linspace(a,b,nx+1)
#c1 = 0.8125
#c2 = 0.511379576154921
#yt = exp(2*xt)*(c1*cos(2*xt)+c2*sin(2*xt)) + 3/8*xt - 5/16
c1 = 13/16
c2 = 15/16 * exp(-2) * 1 / sin(2) - 13/16 * 1 / tan(2)
yt = exp(2*xt)*(c1*cos(2*xt)+c2*sin(2*xt)) + 3/8*xt - 5/16
#print(xt,yt)
plot(xt,yt,'--', linewidth = 2, color = 'black', label = 'Аналитическое решение') 
legend()
show()