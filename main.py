from scipy.interpolate import CubicSpline
from tkinter import *
import math
from math import  cos, sin

stat = 0
dt = 0.001
G = 6.67*10**(-7)
c = 65
r = 0.05
n = 6
r3max = 1e-8

th = 2*r/c
j = int(math.ceil(th/dt))
h =max(j,50) #шаги
km = 350
xc = 600
yc = 300
t = 0


wind = Tk()
layer = Canvas(wind, width=2*xc, height=2*yc)
layer.pack()
inp_dt = StringVar()
inp_n = StringVar()
inp_c = StringVar()
def changeN():
    global n, pr_n
    n = int(inp_n.get())
    layer.delete(pr_n)
    pr_n = layer.create_text(50, 75, text='n = ' + str(n))
def changeC():
    global c, pr_c
    c = int(inp_c.get())
    layer.delete(pr_c)
    pr_c = layer.create_text(50, 100, text='c = ' + str(c))
def changeDt():
    global dt,pr_dt
    dt = float(inp_dt.get())
    layer.delete(pr_dt)
    pr_dt = layer.create_text(50, 25, text='dt = ' + str(dt))
def click_start():
    global stat
    stat = 1
    start()
def click_stop():
    global stat
    stat = 0
def click1():
    global km, pr_km
    if km > 51:km -= 50
    layer.delete(pr_km)
    pr_km = layer.create_text(50, 50, text='Km = ' + str(km))
def click2():
    global km, pr_km
    km += 50
    layer.delete(pr_km)
    pr_km = layer.create_text(50, 50, text='Km = ' + str(km))
def click3():
    global dt, pr_dt
    dt /= 2
    layer.delete(pr_dt)
    pr_dt = layer.create_text(50, 25, text='dt = ' + str(dt))
def click4():
    global dt, pr_dt
    dt *= 2
    layer.delete(pr_dt)
    pr_dt = layer.create_text(50, 25, text='dt = ' + str(dt))
btn1 = Button(wind, text=' km - 50 ', command=click1)
btn2 = Button(wind, text=' km + 50 ', command=click2)
btn3 = Button(wind, text=' dt / 2 ', command=click3)
btn4 = Button(wind, text=' dt x 2 ', command=click4)
btnStart = Button(wind, text=' Start ', command=click_start)
btnStop = Button(wind, text=' Stop ', command=click_stop)
btnDt= Button(text="Change dt",command=changeDt)
btnN = Button(text="Change N",command=changeN)
btnC = Button(text="Change C",command=changeC)
mDt = Entry(textvariable=inp_dt)
mN = Entry(textvariable=inp_n)
mC = Entry(textvariable=inp_c)
btn1.pack(side=LEFT)
btn2.pack(side=LEFT)
btn3.pack(side=LEFT)
btn4.pack(side=LEFT)
btnStart.pack(side=BOTTOM)
btnStop.pack(side=BOTTOM)
mDt.place(relx=.87, rely=.80)
btnDt.place(relx=.87, rely=.83)
mN.place(relx=.87, rely=.90)
btnN.place(relx=.87, rely=.93)
mC.place(relx=.87, rely=.70)
btnC.place(relx=.87, rely=.73)

class line:
    def __init__(self,t0,t1,x_0,x_1):
        self.t0 = t0
        self.x0 = x_0
        self.k = (x_1-x_0)/(t1-t0)
    def __call__(self, t):
        return self.x0+self.k*(t-self.t0)

#Создание класса объекта материальной точки со всеми свойствами

class dot:
    def __init__(self):
        self.ball = layer.create_oval( 3, yc - 6, 15, yc + 6, fill='blue', width=0)
        self.m = 10000.0
        self.x = []
        self.y = []
        self.seg_x = []
        self.seg_y = []
    def __del__(self):
        print("deleted")

#Поиск момента времени, в который точка 1 видит точку 2 с учетом ограничения скорости гравитации методом двоичного поиска
def binSearch(numDot1, numDot2, t):
    #t = round(t,4)
    count = 0
    delt=t
    l = -1
    i = int(t/dt)
    while round(c*delt,2) != round(l,2):
        tau = t - delt
        j = int(tau/dt)
        sgx = dots[numDot2].seg_x[j](tau)
        sgy = dots[numDot2].seg_y[j](tau)
        x02=(dots[numDot1].x[i]- sgx)**2
        y02=(dots[numDot1].y[i]- sgy)**2
        l = math.sqrt(x02 + y02)
        if c*delt > l:
            delt -= delt / 2
        else:
            delt += delt / 2
       # print('l ',l)
    #print(count)
    return delt
#Поиск момента времени, в который точка 1 видит точку 2 с учетом ограничения скорости гравитации методом хорд
def timeSearch(numDot1, numDot2, t):
    tleft = 0
    tright = t - 0.0001*dt
    i = int(t/dt)
    xi = dots[numDot1].x[i]
    yi = dots[numDot1].y[i]
    def compute_l(tau):
        j = int(tau/dt)
        #print(j,len(dots[numDot2].seg_x)
        dx = dots[numDot2].seg_x[j](tau) - xi
        dy = dots[numDot2].seg_y[j](tau) - yi
        return c*(t - tau) - math.sqrt(dx*dx + dy*dy)
    fleft = compute_l(tleft)
    fright = compute_l(tright)
    #count = 0
    tmid = 0.5*(tleft + tright)
    while math.fabs(tleft - tright) > 1e-6:
        tmid = (fleft*tright - fright*tleft) / (fleft - fright)
        if (tmid - tleft ) < 1e-7 or (tright - tmid) < 1e-7:
            tmid = 0.5*(tleft + tright)
        fmid = compute_l(tmid)
        if fleft*fmid < 0:
            tright = tmid
            fright = fmid
        else:
            tleft = tmid
            fleft = fmid
        #count += 1
    #print(count)
    #print(tmid, t)
    return tmid
#Находим ускорение точки в момент времени t
def acc(n1,t):
    a_x = 0
    a_y = 0
    j = int(t/dt)
    for dotn in range(n):
        if dotn != n1:
            t0 = timeSearch(n1,dotn,t)
            b = int(t0/dt)
            #print(t0)
            r3 = (math.sqrt((dots[dotn].seg_x[b](t0) - dots[n1].x[j])**2+
                            (dots[dotn].seg_y[b](t0) - dots[n1].y[j])**2)**3)
            #print(r3,n1,dotn, t-t0)
            r3 = max(r3,r3max)
            a_x += (dots[dotn].m * (dots[dotn].seg_x[b](t0) - dots[n1].x[j]))/r3
            a_y += (dots[dotn].m * (dots[dotn].seg_y[b](t0) - dots[n1].y[j]))/r3

    return G * a_x,G * a_y

#Находим координаты точки в следующий момент времени используя метод Верле
def VerletInt(numDot,t): #xn+1 = 2xn - xn-1 +an*dt^2
    i = int(t/dt)
    Wx,Wy = acc(numDot,t)
    x_t1 = 2*dots[numDot].x[i] - dots[numDot].x[i-1] + Wx*dt**2
    y_t1 = 2*dots[numDot].y[i] - dots[numDot].y[i-1] + Wy*dt**2
    return x_t1,y_t1

#Сплайны для сегмента
def sx(numDot,a):
    dotNum = max(a,4)
    t_s = dots[numDot].time[-4:]
    #print(t_s)
    x_s = [dots[numDot].x[dotNum-4],dots[numDot].x[dotNum-3],dots[numDot].x[dotNum-2],dots[numDot].x[dotNum-1]]
    x_tck = CubicSpline(t_s,x_s)
    #print(x_tck(0))
    return x_tck
def sy(numDot,a):
    dotNum = max(a,4)
    t_s = dots[numDot].time[-4:]
    y_s = [dots[numDot].y[dotNum-4],dots[numDot].y[dotNum-3],dots[numDot].y[dotNum-2],dots[numDot].y[dotNum-1]]
    y_tck = CubicSpline(t_s,y_s)
    #print(y_tck(0))
    return y_tck
def velocity(x0,x1):
    return  (x1 - x0)/dt

dots=[]
#Создаем материальные точки
def start():
    global dt,G,c,r,n,r3max,th,j,h,km,xc,yc,t,pr_km,pr_dt,pr_n,layer,pr_c
    layer.delete('all')
    layer.create_line(xc, 0, xc, 2*yc, width=1)
    layer.create_line(0, yc, 2*xc, yc, width=1)
    layer.create_oval( xc-3, yc - 3, xc+3, yc + 3, fill='red', width=0)

    for a in range(n):
        dot_new = dot()
        dots.append(dot_new)
        for b in range(h):
            t = b * dt
            arc = t + a
            dot_new.x.append(r*cos(arc))
            dot_new.y.append(r*sin(arc))
            if b > 0:
                dot_new.seg_y.append(line(t-dt, t, dot_new.y[b-1], dot_new.y[b]))
                dot_new.seg_x.append(line(t-dt, t, dot_new.x[b-1], dot_new.x[b]))
    pr_km = layer.create_text(50,50,text='Km = ' + str(km))
    pr_dt = layer.create_text(50,25,text='dt = ' + str(dt))
    pr_n  = layer.create_text(50, 75, text='n = ' + str(n))
    pr_c  = layer.create_text(50, 100, text='c = ' + str(c))
    pr_maxV = layer.create_text(70, 125, text='max V = ' + str(0))

    for i in range (h-1,50000): #i шаг, t текущее время
        t = dt * i
        for j in range(n):
            if stat == 0:
                for k in range(len(dots)): del dots[0]
                return 0
            x1,y1=VerletInt(j,t)
            dots[j].x.append(x1)
            dots[j].y.append(y1)
            dots[j].seg_x.append(line(t,t+dt,dots[j].x[i],x1))
            dots[j].seg_y.append(line(t,t+dt,dots[j].y[i],y1))
            layer.moveto(dots[j].ball, 5*x1*km-7+xc,5*y1*km+yc-7)
            layer.create_oval(5*x1*km + xc, 5*y1*km + yc, 5*x1*km + xc, 5*y1*km + yc, width=1)  # траектория
        layer.delete(pr_maxV)
        maxV = max(math.hypot(velocity(dots[j].x[i - 1], dots[j].x[i]), velocity(dots[j].y[i - 1], dots[j].y[i])) for j in range(n))
        pr_maxV = layer.create_text(60, 125, text='max V = ' + str(round(maxV,7)))
        wind.update()
wind.mainloop()
