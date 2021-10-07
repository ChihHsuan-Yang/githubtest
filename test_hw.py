import numpy as np
import math

def bisection(f,a,b,TOL,MaxIters):
    fa = f(a)
    fb = f(b)
    assert fa*fb < 0
    assert b > a
    n = 0
    while n < MaxIters:
        n += 1
        alpha_n = 0.5*(a + b)
        f_alpha_n = f(alpha_n)
        if 0.5*(b - a) < TOL and np.abs(f_alpha_n) < TOL:
            break
        else:
            if fa * f_alpha_n < 0:
                b = alpha_n
            else:
                a = alpha_n        
    return alpha_n,n,f_alpha_n

def falsepos(f,a,b,TOL,MaxIters):
    alpha = []
    fa = f(a)
    fb = f(b)
    
    assert fa*fb < 0
    assert b > a
    
    n = 1
    alpha_1 = a - fa*(b-a)/(fb - fa)
    f_alpha_1 = f(alpha_1)
    alpha.append(alpha_1)
    
    if f_alpha_1*fa < 0:
        b = alpha_1
        fb = f_alpha_1
    else:
        a = alpha_1
        fa = f_alpha_1
    
    while n < MaxIters:
        n += 1
        alpha_n = a - fa*(b-a)/(fb - fa)
        alpha.append(alpha_n)
        
        f_alpha_n = f(alpha_n)
        
        if f_alpha_n*fa < 0:
            b = alpha_n
            fb = f_alpha_n
        else:
            a = alpha_n
            fa = f_alpha_n
        
        if np.abs(alpha[-1] - alpha[-2]) < TOL and np.abs(f_alpha_n) < TOL:
            break
               
    return alpha_n,n,f_alpha_n


def newton(f,df,alpha_0,TOL,MaxIters):
    n = 0
    alpha = []
    alpha.append(alpha_0)
    fval = f(alpha_0)
    fpval = df(alpha_0)
    
    while n < MaxIters:
        n += 1
        delta = - fval/fpval
        alpha_n = alpha[-1] + delta
        alpha.append(alpha_n)
        
        fval = f(alpha_n)
        
        if np.abs(delta) < TOL and np.abs(fval) < TOL:
            break
        else:
            fpval = df(alpha_n)
        
    return alpha_n,n,fval

def secant(f,alpha_0,alpha_1,TOL,MaxIters):
    alpha = []
    alpha.append(alpha_0)
    alpha.append(alpha_1)
    
    n = 1
    f_alpha_n_2 = f(alpha_0)
    f_alpha_n_1 = f(alpha_1)
    while n < MaxIters:
        n += 1
        alpha_n = alpha[-1] -  f_alpha_n_1*((alpha[-1] - alpha[-2])/(f_alpha_n_1 - f_alpha_n_2))
        alpha.append(alpha_n)
        f_alpha_n = f(alpha_n)
        if np.abs(alpha[-1] - alpha[-2]) < TOL and np.abs(f_alpha_n) < TOL:
            break
        else:
            f_alpha_n_2 = f_alpha_n_1
            f_alpha_n_1 = f_alpha_n
    
    return alpha_n,n,f_alpha_n

def f1(x):
    y=np.sin(x)-x-1
    return y
def f2(x):
    y=x*(1-np.cos(x))
    return y
def f3(x):
    y=np.exp(x)-x**2+3*x-2
    return y
def df1(x):
    y=np.cos(x)-1
    return y
def df2(x):
    y=(1-np.cos(x))+x*np.sin(x)
    return y
def df3(x):
    y=np.exp(x)-2*x+3
    return y
g1 = lambda x: (np.sin(x) - x - 1)**2

TOL = 1e-12
MaxIters = 500
a = -2
b = 1
x_0 = 1
x_1 = 0.9

print('for equation 1')
print('bisection : ')
print(bisection(f1,a,b,TOL,MaxIters))
print('falsepos : ')
print(falsepos(f1,a,b,TOL,MaxIters))
print('newton : ')
print(newton(f1,df1,x_0,TOL,MaxIters))
print('secant : ')
print(secant(f1,x_0,x_1,TOL,MaxIters))
print('\n')

print('for equation 2')
print('bisection : ')
print(bisection(f2,a,b,TOL,MaxIters))
print('falsepos : ')
print(falsepos(f2,a,b,TOL,MaxIters))
print('newton : ')
print(newton(f2,df2,x_0,TOL,MaxIters))
print('secant : ')
print(secant(f2,x_0,x_1,TOL,MaxIters))
print('\n')

print('for equation 3')
print('bisection : ')
print(bisection(f3,a,b,TOL,MaxIters))
print('falsepos : ')
print(falsepos(f3,a,b,TOL,MaxIters))
print('newton : ')
print(newton(f3,df3,x_0,TOL,MaxIters))
print('secant : ')
print(secant(f3,x_0,x_1,TOL,MaxIters))
print('\n')


# The question 6:
R = 0.08205
T = 313
P = 2 
a = 12.87
b = 0.1142
n = 1

equation = lambda V : (P + a/V**2)*(V - b) - n*R*T

minV = 0.1
maxV = 100
TOL = 1e-12
MaxIters = 500
bisection(equation,0.1,100,TOL,MaxIters)

# The question 7:
r = 5
p = 1
z = 0.6 

equation = lambda h : 1/3 *np.pi* (3*r*h**2 -h**3)*p -4/3 * np.pi*r**3 *z

minV = -100
maxV = 100
TOL = 1e-12
MaxIters = 500

numInterval = int((maxV-minV)*2+1)
intervalRange = np.linspace(minV ,maxV,numInterval)

fa = equation(intervalRange[0])
fb = equation(intervalRange[1])

for i in range(1,numInterval-1 ):
    if fa*fb < 0:
        alpha_,n,f_alpha_ =  falsepos(equation,intervalRange[i-1],intervalRange[i],TOL,MaxIters)
        print('Find a root : ' + str(alpha_))
        print('The function value at the root : ' + str(f_alpha_))

    fa = fb
    fb = equation(intervalRange[i+1])
        
# problem 5 c   
print('for discussion (c) :')
print(bisection(g1,a,b,TOL,MaxIters))  
