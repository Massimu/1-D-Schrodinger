#Computational physics- written by tao pang - Chapter 4 - section 9
#An example of solving the eigenvalue problem of the one-dimensional Schroedinger equation via the secant and Numerov methods..

#packages
import matplotlib.pyplot as plt
import numpy as np

#initial values for definig equation and calculation
nx=500
m=10
ni=10
delta=1e-6
e=2.4
de=0.1
x1=-10
x2=10
h=float((x2-x1)/500.)
#initial values 
ql=[{}]*(nx+1)
qr=[{}]*(nx+1)
s=[0]*(nx+1)
u=[{}]*(nx+1)
#------------------------------Simpson method for finding the estimated amount for our initial values--------------------------
def simpson(y,h):
    n=int(len(y)-1)
    s0=0
    s1=0
    s2=0
    for i in range(1,n):
        s0 += y[i]
        s1 += y[i-1]
        s2 += y[i+1]
    s=float((s1+4*s0+s2)/3.)
    if (n+1)%2==0:
        return h*(s+(5*y[n]+8*y[n-1]-float(y[n-2]/12.)))
    else:
        return h*s
#----------------------------------------------------------------
#Method to provide the given potential in the problem.
#----------------------------Potential---------------------------
import math 
def v(x):
    alpha = 1
    lamb =4
    return alpha*alpha*lamb*(lamb-1)*float(0.5-float(1/(math.pow(math.cosh(alpha*x),2))))/2.
#---------------------------------------------------------------
#----------------------------Secant method : Find the eigenvalue via the secant search-----------------------------
def Sacant(n,delta,x,dx):        
    k=0
    x1=x+dx
    while abs(dx)>delta and k<n :
           
            x2 = x1-f(x1)*float(((x1-x)/(f(x1)-f(x))))
            x=x1
            x1=x2
            dx=x1-x
            k+=1
    return x1
#-------------------------------------------------------------------
#-------------------------------Numerov method-----------------------------
#which is an extremely accurate scheme for linear differential equations without the first-order derivative term
#defined in textbook
def numerov(m,h,u0,u1,q,s):
    u=[{}]*m
    u[0]=u0
    u[1]=u1
    g=float((h*h)/12)
    for i in range(1,m-1):
        c0 = 1+g*q[i-1]
        c1 = 2-10*g*q[i]
        c2 = 1+g*q[i+1]
        d = g*(s[i+1]+s[i-1]+10*s[i])
        u[i+1] = float((c1*u[i]-c0*u[i-1]+d)/c2)  
    return u
#-------------------------------------------------------------------          
#definition the function which is suppose to calculate its root to find the matching points between potential and wavefunction
#----------------------------------F(x)-----------------------------
def f(E):
    y=[{}]*(nx+1)
    u0=0
    u1=0.01
    for i in range(0,nx+1):
        x = -10+i*h
        ql[i] = 2*(E-v(x))
        qr[nx-i] = ql[i]
   
    im=0
    for i in range(nx):
        if (ql[i]*ql[i+1]<0 and ql[i]>0): 
            im=i
        
    nl=im+2
    nr=nx-im+2
    c0 = 0
    c1 = 0
    c2 = 0 
    ul=[{}]*nl
    ul[0]=0
    ul[1]=0.01
    g=float((h*h)/12)
    for i in range(1,nl-1):
        c0 = 1+g*ql[i-1]
        c1 = 2-10*g*ql[i]
        c2 = 1+g*ql[i+1]
        d = g*(s[i+1]+s[i-1]+10*s[i])
        ul[i+1] =float((c1*ul[i]-c0*ul[i-1]+d)/c2)
    
    
    c0 = 0
    c1 = 0
    c2 = 0    
    ur=[{}]*nr
    ur[0]=u0
    ur[1]=u1
    g=float((h*h)/12)
    for i in range(1,nr-1):
        c0 = 1+g*qr[i-1]
        c1 = 2-10*g*qr[i]
        c2 = 1+g*qr[i+1]
        d = g*(s[i+1]+s[i-1]+10*s[i])
        ur[i+1] = float((c1*ur[i]-c0*ur[i-1]+d)/c2)  

        
    
    #left Wave function
    ratio = float(ur[nr-2]/ul[im])
    #print ratio
    for i in range(im+1):
        u[i] = ratio*ul[i]
        y[i] = u[i]*u[i]
        
    #right wave function
    for i in range(nr-1):
        u[i+im] = ur[nr-i-2]
        y[i+im] = u[i+im]*u[i+im]
        
    #normalize the wave function
    
    Sum=simpson(y,h)
    Sum = math.sqrt(Sum)
    #print Sum
    for i in range(nx+1):
        u[i] /= Sum 
    f0 = ur[nr-1]+ul[nl-1]-ur[nr-3]-ul[nl-3]
    return float(f0/(2*h*ur[nr-2])) 
    
#main section
e = Sacant(ni,delta,e,de)

#plot potetial and wave function :
x=x1
V = [{}]*(nx+1)
alpha = 1
lamb =4
for i in range(nx+1):
   x    = x1 + i*h
   V[i] =   alpha*alpha*lamb*(lamb-1)*float(0.5-float(1/(math.pow(math.cosh(alpha*x),2))))/2.
    




x=np.linspace(-10,10,501)
plt.plot(x,u,x,V)
plt.grid()
plt.title("The eigenvalue Schrodinger equation")
plt.xlabel("$x$")
plt.ylabel("$V(x),U(x)$")
plt.show()

    

        
    


 
