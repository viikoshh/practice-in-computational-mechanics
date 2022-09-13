from numpy import array, linspace, zeros, sin, trapz, cos, exp, tan
from numpy.linalg import solve
from matplotlib.pylab import figure, plot, show, legend

def integ(a,b,m):
    x = linspace(a,b)
    if m: l = x - a
    else: l = b - x
    return trapz(l/(b-a)*(3 * x - 4), x)
    
def mke(n,a,b,ya,yb):
    h = (b - a) / n  #; print(h)
    E = n + 1 #количество узлов
    A = zeros((E,E))
    B = zeros((E))
    x = [i*h for i in range(n+1)] #;print(x)
    Ak =  1 / h *array([[1, -1], [-1, 1]]) + 2 * array([[-1, 1], [-1, 1]])\
            - 8 * h * array([[1/3, 1/6], [1/6, 1/3]])
    #print(Ak)
    for i in range(E-1):
        A[i:i+2,i:i+2] += Ak
    #print(A)
    for i in range(E-1):
        x_k = x[i]; x_k1 = x[i+1]
        Bk = [0, 0]
        Bk[0] = integ(x_k, x_k1, 0)
        Bk[1] = integ(x_k, x_k1, 1)
        B[i:i+2] += Bk
    #print(B)
    Aa = A[1:n,1:n] #;print(Aa)
    Bb = -B[1:n]  #;print(Bb) 
    Bb[0] -= A[1,0]*ya; Bb[n-2] -= A[n-1,n]*yb
    #print(Aa); print(Bb)
    ab = solve(Aa,Bb)  #; print(ab)
    
    X = zeros((n+1));X[0] = a;X[n] = b
    Y = zeros((n+1));Y[0] = ya;Y[n] = yb
    for i in range(1,n):
        X[i] = x[i]
        Y[i] = ab[i-1]
    #print(X);print(Y)
    return X,Y

sn = (3, 6, 100) 
a, b = 0, 1
ya, yb = 0.5, 1
rgb = ('r','g','b')
#line = ('-','--',':')
#linew = (1, 2, 4) 
figure()
for k in range(len(sn)):
    n = sn[k]
    x, y = mke(n,a,b,ya,yb)
    #plot(x,y,linestyle = line[k],linewidth = linew[k], color = rgb[k], label = f'N = {n}')
    plot(x,y, color = rgb[k], label = f'N = {n}')
#точное решение
nt = 1000
xt = linspace(a,b,nt+1)
c1 = 13/16
c2 = 15/16 * exp(-2) * 1 / sin(2) - 13/16 * 1 / tan(2)
yt = exp(2*xt)*(c1*cos(2*xt)+c2*sin(2*xt)) + 3/8*xt - 5/16
#print(xt,yt)
plot(xt,yt,'--', linewidth = 2, color = 'black', label = 'Аналитическое решение')     
legend()
show() 