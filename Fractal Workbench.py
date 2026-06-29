import Fractal_IO as io
import Fractal_Tools as tools

import time as time
import datetime
import math
import matplotlib.pyplot as plt

def print_date_and_time():
    t = time.localtime()
    d = datetime.date.today()
    print(time.strftime("%H:%M:%S", t), d)
    
def process_image():
    filename = "image15E.tif"
    x_min, x_max = 0, 542
    y_min, y_max = 0, 542
    
    lo, hi = .11, 1
    xs, ys, ms = tools.mass_fractalizer(filename, lo, hi, x_min, x_max, y_min, y_max, trunc = 4)
    write_filename = filename[0:-4]+"_"+str(lo)+"_"+str(hi)+".txt"
    #write_filename = "tri_mass.txt"
    io.write_mass_fractal("sprite1_11p.txt", xs, ys, ms)
    
    #filename = "sprite1.tif"
    #x_min, x_max = 0, 0
    #y_min, y_max = 0, 0
    
    lo, hi = .12, 1
    xs, ys, ms = tools.mass_fractalizer(filename, lo, hi, x_min, x_max, y_min, y_max, trunc = 4)
    write_filename = filename[0:-4]+"_"+str(lo)+"_"+str(hi)+".txt"
    #write_filename = "tri_mass.txt"
    io.write_mass_fractal("sprite1_12p.txt", xs, ys, ms)
    
    filename = "image21E.tif"
    x_min, x_max = 0, 652
    y_min, y_max = 0, 652
    
    lo, hi = .14, 1
    xs, ys, ms = tools.mass_fractalizer(filename, lo, hi, x_min, x_max, y_min, y_max, trunc = 4)
    write_filename = filename[0:-4]+"_"+str(lo)+"_"+str(hi)+".txt"
    #write_filename = "tri_mass.txt"
    io.write_mass_fractal("sprite2_14p.txt", xs, ys, ms)
    
    #filename = "image21E.tif"
    #x_min, x_max = 0, 652
    #y_min, y_max = 0, 652
    
    lo, hi = .16, 1
    xs, ys, ms = tools.mass_fractalizer(filename, lo, hi, x_min, x_max, y_min, y_max, trunc = 4)
    write_filename = filename[0:-4]+"_"+str(lo)+"_"+str(hi)+".txt"
    #write_filename = "tri_mass.txt"
    io.write_mass_fractal("sprite2_16p.txt", xs, ys, ms)
    
def process_points():
    
    #filenames = ["sprite1_11p.txt", "sprite1_12p.txt", "sprite2_14p.txt", "sprite2_16p.txt"]
    filenames = ["sprite2_14p.txt", "sprite2_16p.txt"]
    #filenames = ["tri_mass.txt"]
    
    Qs = [Q/2 for Q in range(-30, 30+1)]
    
    for filename in filenames:
        xs, ys, ms = io.read_mass_fractal(filename)
        zs = [0 for x in xs]
        
        num_rs = 100
        max_r = math.sqrt((max(xs)-min(xs))**2 + (max(ys)-min(ys))**2)
        min_r = 4
        Rs = [min_r*math.pow((max_r/min_r),(i/num_rs)) for i in range(num_rs+1)]
        
        print_date_and_time()
        
        lnCs = tools.Correlation_Integrator(xs, ys, zs, ms, Qs, Rs)
        lnRs = [math.log(R) for R in Rs]
        
        write_filename = filename[0:-4] + "_lnC(q,r).txt"
        io.write_XY2file(write_filename, Qs, lnRs, lnCs)
        
    print_date_and_time()
        

def fitting():
    gen_log = True
    swap = -1

    lower, upper = 2, 5
    Qs, xs, yss = io.read_XY("sprite1_11p_lnC(q,r).txt")
    ms5, paramss = tools.L_sloper(xs, yss, Qs, lower, upper, gen = gen_log, swap = swap)
    io.plot_L_fit_lines(Qs, ms5, paramss, xs, yss, upper, lower, title = 'A')
    
    Qs, xs, yss = io.read_XY("sprite1_12p_lnC(q,r).txt")
    ms6, paramss = tools.L_sloper(xs, yss, Qs, lower, upper, gen = gen_log, swap = swap)
    io.plot_L_fit_lines(Qs, ms6, paramss, xs, yss, upper, lower, title = 'C')
    
    Qs, xs, yss = io.read_XY("sprite2_14p_lnC(q,r).txt")
    ms7, paramss = tools.L_sloper(xs, yss, Qs, lower, upper, gen = gen_log, swap = swap)
    io.plot_L_fit_lines(Qs, ms7, paramss, xs, yss, upper, lower, title = 'B')
    
    Qs, xs, yss = io.read_XY("sprite2_16p_lnC(q,r).txt")
    ms8, paramss = tools.L_sloper(xs, yss, Qs, lower, upper, gen = gen_log, swap = swap)  
    io.plot_L_fit_lines(Qs, ms8, paramss, xs, yss, upper, lower, title = 'D')
    
    Qs, xs, yss = io.read_XY("HD_tri_mass_lnC(q,r).txt")
    lo, hi = 3, 5
    ms9, bs = tools.sloper(xs, yss, lower = lo, upper = hi, fit_C_range = False)
    ms10 = tools.analytical_multi_sierp_spec(Qs, ps = (.20, .38, .42))
    
    plt.figure()
    
    color = "green"
    plt.fill_between(Qs, ms5, ms6, alpha = .5, color = color, label = "Sprite 1")
    plt.plot(Qs, ms5, color = color)
    plt.plot(Qs, ms6, color = color)
    
    color = "blue"
    plt.fill_between(Qs, ms7, ms8, alpha = .5, color = color, label = "Sprite 2")
    plt.plot(Qs, ms7, color = color)
    plt.plot(Qs, ms8, color = color)
    
    
    plt.plot(Qs, ms9, alpha = .5, color = 'red', label = "Multifractal Sierpinski Triangle (Numerical)")
    plt.plot(Qs, ms10, alpha = .5, color = 'orange', label = "Multifractal Sierpinski Triangle (Analytical)")
    
    
    plt.legend()
    plt.xlabel("q")
    plt.ylabel("D(q)")
    plt.gca().set_ylim(bottom = 0)
    plt.tight_layout()
    
    
def main():
    #This is for playing around with thresholds and generating the data
    #process_image()
    #This calculates C_q(r) (very slow)
    #process_points()
    #This actually generates the spectrum (quick)
    fitting()

main()