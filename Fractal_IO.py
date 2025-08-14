import math
import matplotlib.pyplot as plt
import numpy as np                                                               

def write_mass_fractal(filename, xs, ys, ms):
    print("Writing:", filename)
    file = open(filename, 'w')
    for i in range(len(xs)):
        point = str(xs[i])+" "+str(ys[i])+" "+str(ms[i])+"\n"
        file.write(point)
    file.close()
    
def read_mass_fractal(filename):
    xs = []
    ys = []
    ms = []
    file = open(filename, 'r')
    line = file.readline()
    while(line != ''):
        x, y, m = line.split()
        xs.append(eval(x))
        ys.append(eval(y))
        ms.append(eval(m))
        line = file.readline()
    return xs, ys, ms

def write_XY2file(filename, Qs, Xs, Yss):
    print("Writing:", filename)
    file = open(filename, 'w')
    
    for q in Qs:
        file.write(str(q)+" ")
    file.write("\n")
    
    for x in Xs:
        file.write(str(x)+" ")
    file.write("\n")
    
    for i in range(len(Yss)):
        Ys = Yss[i]
        for y in Ys:
            file.write(str(y)+" ")
        file.write("\n")
        
    file.close()
    
def read_XY(filename):
    qs, xs, yss = [], [], []
    file = open(filename, 'r')
    
    line = file.readline()
    line_eles = line.split()
    for q in line_eles:
        qs.append(eval(q))
        
    line = file.readline()
    line_eles = line.split()
    for x in line_eles:
        xs.append(eval(x))
        
    line = file.readline()
    while(line != ''):
        ys = []
        line_eles = line.split()
        for y in line_eles:
            ys.append(eval(y))
        yss.append(ys)
        line = file.readline()
    return qs, xs, yss

def plot_fit_lines(Qs, Ds, bs, Xs, Yss, lower = 0, upper = 0, title = "Generic Fit Lines", plot_all_points = True, fit_C_range = True):
    valid_Qs = [-15, -10, -5, 0, 5, 10, 15]
    #valid_Qs = [-13, -12, -11, -10, -9]
    #plot the data
    plt.figure(title)
    if(not fit_C_range):
        lower_ind, upper_ind = len(Xs)-1, 0
        for i in range(len(Xs)):
            x = Xs[i]
            if(x > lower and x < upper):
                if i < lower_ind:
                    lower_ind = i
                if i > upper_ind:
                    upper_ind = i
        print(lower_ind, upper_ind)
        print(Xs[lower_ind], Xs[upper_ind])
    for i in range(len(Qs)):
        Q = Qs[i]
        if Q in valid_Qs:
            Ys = Yss[i]
            lower_ind, upper_ind = len(Ys)-1, 0
            for j in range(len(Ys)):
                if(fit_C_range):#then we are fitting to a correlation range
                    c = math.exp(Ys[j])#this isn't true, is it? it is.
                    if(c >= lower and c <+ upper):# if c is in the range
                        if(j < lower_ind):
                            lower_ind = j #this will only trigger once
                        if(j > upper_ind):
                            upper_ind = j #This will trigger so many times
            #the data
            if(plot_all_points):
                line = plt.plot(Xs, Ys, ms = 1.5, marker = 'd',  linestyle = "none")
                c = line[0].get_color()
                plt.plot([],[], marker = 'd', color = c, linestyle = "none", label = str(Q))
            else:
                line = plt.plot(Xs[lower_ind:upper_ind], Ys[lower_ind:upper_ind], ms = 1.5, marker = 'd',  linestyle = "none")
                c = line[0].get_color()
                plt.plot([],[], marker = 'd', color = c, linestyle = "none", label = str(Q))
            #the fit lines
            m, b = Ds[i], bs[i]
            if(fit_C_range):
                x0 = (Ys[lower_ind]-b)/m
                xf = (Ys[upper_ind]-b)/m
                #x0 = (math.log(lower)-b)/m
                #xf = (math.log(upper)-b)/m
            else:
                #x0 = Xs[lower_ind]
                #xf = Xs[upper_ind]
                x0 = lower
                xf = upper
            plt.plot([x0, xf], [m*x0+b, m*xf+b], color = c, linestyle = "solid")
    plt.legend()
    plt.xlabel("ln(r)")
    plt.ylabel("ln(C(q,r)")
    plt.title(title)
    
def logistic(xs, L, k, a, b):
    return (L-b)/(1+np.exp(k*(xs-a)))+b

def gen_logistic(xs, L, k, a, b, v):
    v = math.exp(v)
    return (L-b)/np.power(1+np.exp(k*(xs-a)), v)+b
    
def plot_L_fit_lines(Qs, Ds, paramss, Xs, Yss, upper = 0, lower = 0, title = "title"):
    valid_Qs = [-15, -10, -5, -1, 0, 1, 5, 10, 15]
    #valid_Qs = [-15, -10, -5, -3, -1]
    #valid_Qs = [-15, 0, 15]
    #valid_Qs = [-15, -11, -7, -3, 0, 3, 7, 11, 15]
    #valid_Qs = [-2, -1.5, -1, -.5, 0, 3, 7, 11, 15]
    if(upper != lower):
        lower_ind, upper_ind = len(Xs)-1, 0
        for i in range(len(Xs)):
            x = Xs[i]
            if(x > lower and x < upper):
                if i < lower_ind:
                    lower_ind = i
                if i > upper_ind:
                    upper_ind = i
    
    #plt.figure(title, figsize = [8, 8])
    #plt.figure(figsize = [3, 12])
    
    fig, axes = plt.subplots(nrows = 3, ncols = 3, sharex = True, sharey = True, figsize = [8,8])
    #fig, axes = plt.subplots(1, 3, sharex = False, sharey = False, num = title)
    counter = 0
    for i in range(len(Qs)):
        if(Qs[i] in valid_Qs):
            #ax = axes[counter]
            ax = axes[counter%3][counter//3]
            counter += 1
            if(len(paramss[i]) == 5):
                L, k, a, b, v = paramss[i]
                ax.plot(Xs, gen_logistic(Xs, L, k, a, b, v), color = "green")
                #b = y0 - mx0
                v = math.exp(v)
                b = (L-b)/math.pow(1 + 1/v, v) + b - (a - math.log(v)/k)*Ds[i]
                y0, yf = min(Yss[i]), max(Yss[i]) 
                x0, xf = (y0-b)/Ds[i], (yf-b)/Ds[i] 
                ax.plot([x0, xf], [y0, yf], color = "orange", label = "q = " + str(int(Qs[i]))) 
            elif(len(paramss[i]) == 4): 
                L, k, a, b = paramss[i] 
                ax.plot(Xs, logistic(Xs, L, k, a, b), color = "green") 
                b = -a*Ds[i] + (L + b)/2 
                y0, yf = min(Yss[i]), max(Yss[i]) 
                x0, xf = (y0-b)/Ds[i], (yf-b)/Ds[i] 
                ax.plot([x0, xf], [y0, yf], color = "orange", label = "q = " + str(int(Qs[i]))) 
            elif(len(paramss[i]) == 1):
                b, = paramss[i]
                if(upper != lower):
                    x0, xf = Xs[lower_ind], Xs[upper_ind]
                else:
                    x0, xf = Xs[0], Xs[-1]
                ax.plot([x0, xf], [Ds[i]*x0+b, Ds[i]*xf+b], color = "orange", label = "q = " + str(int(Qs[i])))
            ax.plot(Xs, Yss[i], ms = 1.5, marker = '|',  linestyle = "none", color = "blue")
            ax.legend()
            
            #ax.set_title("q = " + str(int(Qs[i])))
            #ax.set_title(str(counter%3) + "," + str(counter//3))
            if((counter-1)%3 == 2):
                ax.set_xlabel("ln(r)", size = 'x-large')#, fontsize = "x-small"
            if((counter-1)//3 == 0):
                ax.set_ylabel("ln(C(q,r)", size = 'x-large')
            
            """
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange((start+.5)//1, (end+.5)//1, 1), fontsize = "xx-small")
            start, end = ax.get_ylim()
            ax.yaxis.set_ticks(np.arange((start+.5)//1, (end+.5)//1, 2), fontsize = "xx-small")
            """
    #fig.suptitle(title)
    plt.suptitle(title, size = 'xx-large')
    plt.tight_layout()
    