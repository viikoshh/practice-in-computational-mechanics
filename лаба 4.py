from numpy import zeros, sin, cos, linspace, exp, tan
from numpy.linalg import solve
from sympy import diff, symbols, lambdify
from matplotlib.pylab import figure, plot, show, legend

def coll(n,a,b,ya,yb,nx):
    h = (b-a)/n #шаг разностной сетки
    xi = [i*h for i in range(1,n)] #;print(xi)
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
            #A[i,j-1] = ddphi(xx, j)-4*(0.5+dphi(xx, j))+8*(0.5*xx+0.5+phi(xx, j))
            A[i,j-1] = ddphi(xx, j) - 4 * dphi(xx, j) + 8 * phi(xx, j)
        #B[i] = 3 * xx - 4
        B[i] =  -xx - 6
    #print(A);print(B);
    ab = solve(A,B)   #;print(ab)
    
    x = linspace(a,b,nx)
    y = zeros(nx);y[0] = ya; y[nx-1] = yb
    for i in range(1,nx-1):
        xx = x[i] #;print(xx)
        y[i] = 0.5*(xx + 1)
        for j in range(1,n):
           y[i] += ab[j-1] * xx ** j * (1 - xx)
    #print(x); print(y)
    return x,y
    
#d^2y/dx^2 = sin(pi*x), y(0)=1, y(1)=2
sn = (3, 4, 20) 
a, b = 0, 1
ya, yb = 0.5, 1
rgb = ('r', 'g', 'b')
#line = ('-', '--', ':')
#linew = (1, 2, 4)
nx = 100
figure()
for k in range(len(sn)):
    n = sn[k]
    x, y = coll(n,a,b,ya,yb,nx)
    #plot(x,y,linestyle = line[k],linewidth = linew[k], color = rgb[k])
    plot(x,y, color = rgb[k], label = f'N = {n}')
#точное решение
nt = 100
xt = linspace(a,b,nt+1)
c1 = 13/16
c2 = 15/16 * exp(-2) * 1 / sin(2) - 13/16 * 1 / tan(2)
yt = exp(2*xt)*(c1*cos(2*xt)+c2*sin(2*xt)) + 3/8*xt - 5/16
#print(xt,yt)
plot(xt,yt,'--', linewidth = 2, color = 'black', label = 'Аналитическое решение') 
legend()
show()
