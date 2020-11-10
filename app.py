# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 12:26:57 2020

@author: santi
"""

from PyQt5.QtWidgets import (QApplication, QMainWindow)
from gui import Ui_MainWindow
import numpy as np
from scipy import signal
import time


app = QApplication([])
MainWindow = QMainWindow()
ui = Ui_MainWindow()
ui.setupUi(MainWindow)
MainWindow.show()

def generar(t, tipo, amp, per, duty, primer, width, fs):
    if tipo == "Sinusoide":
        return amp * np.cos(2*np.pi*t/per)
    if tipo == "Cuadrada":
        return amp * signal.square(2*np.pi*t/per, duty/100)
    if tipo == "Sawtooth":
        return amp * signal.sawtooth(2*np.pi*t/per, width)
    if tipo == "Sinc":
        return amp * np.sinc(t/primer)
    if tipo == "Pulso rectangular":
        y = np.zeros(len(t))        
        inicio = int((0-t[0]*fs)) - int(width*fs/2) 
        fin = int((0-t[0]*fs)) + int(width*fs/2)
        y[inicio:fin] = 1
        return amp * y     
    if tipo == "Pulso triangular":
        y = np.zeros(len(t))        
        inicio = int((0-t[0]*fs)) - int(width*fs/2) 
        fin = int((0-t[0]*fs)) + int(width*fs/2)
        if (fin - inicio) % 2 == 0:
            ramp = np.linspace(0,1,(fin - inicio)//2)
            tr = np.hstack([ramp, ramp[::-1]])                             
        else:
            ramp = np.linspace(0,1,(fin - inicio)//2)
            tr = np.hstack([ramp, ramp[::-1]])
        y[inicio:fin] = tr
        return amp * y
    if tipo == "Escalon":
        return np.heaviside(t,1)
    if tipo == "Signo":
        return np.sign(t)
    if tipo == "Exponencial decreciente":
        return np.exp(-width*t) * np.heaviside(t,1)
    if tipo == "Exponencial bilateral":
        return np.exp(-width*abs(t))

    
    
def manipulardg(t0, fs, giro, x):
    y = x
    desp = np.zeros(int(abs(t0*fs)))
    if t0 > 0:
        y = np.concatenate((desp, y), axis=None)
        y = y[:len(x)]
    if t0 < 0:
        y = np.concatenate((y, desp), axis=None)
        y = y[len(desp):]
    if giro:
        y = y[::-1]
        
    return y

def manipulargd(t0, fs, giro, x):
    y = x
    desp = np.zeros(int(abs(t0*fs)))
    if giro:
        y = y[::-1]
    if t0 > 0:
        y = np.concatenate((desp, y), axis=None)
        y = y[:len(x)]
    if t0 < 0:
        y = np.concatenate((y, desp), axis=None)
        y = y[len(desp):]

    return y
        
        
    


def procesar():
    fs = float(ui.fs.text())
    ti = float(ui.t_i.text())
    ts = float(ui.t_s.text())
    T = 1/fs
    t = np.arange(ti, ts + T, T)
    # variables señal 1
    tipo1 = ui.s1_type.currentText()
    amp1 = float(ui.amp_s1.text())
    per1 = float(ui.per_s1.text())
    duty1 = float(ui.duty_s1.text())
    primer1 = float(ui.primer_s1.text())
    width1 = float(ui.width_s1.text())
    # variables señal 2
    tipo2 = ui.s2_type.currentText()
    amp2 = float(ui.amp_s2.text())
    per2 = float(ui.per_s2.text())
    duty2 = float(ui.duty_s2.text())
    primer2 = float(ui.primer_s2.text())
    width2 = float(ui.width_s2.text())
    #genero y grafico señal 1
    y1 = generar(t, tipo1, amp1, per1, duty1, primer1, width1, fs)
    ui.mpl_s1.canvas.axes.clear()
    ui.mpl_s1.canvas.axes.set_title('Se\u00F1al 1 Original')
    ui.mpl_s1.canvas.axes.plot(t, y1)
    ui.mpl_s1.canvas.draw() 
    #genero y grafico señal 2
    y2 = generar(t, tipo2, amp2, per2, duty2, primer2, width2, fs)
    ui.mpl_s2.canvas.axes.clear()
    ui.mpl_s2.canvas.axes.set_title('Se\u00F1al 2 Original')
    ui.mpl_s2.canvas.axes.plot(t, y2)
    ui.mpl_s2.canvas.draw()
    #variables manipular
    giro = ui.giro.isChecked()
    desp_giro = ui.desp_giro.isChecked()
    t0 = float(ui.desp.text())
    #genero y grafico la señal 1 manipulada
    if desp_giro:
        y1m = manipulardg(t0, fs, giro, y1)
    else:
        y1m = manipulargd(t0, fs, giro, y1)
    ui.mpl_s1m.canvas.axes.clear()
    ui.mpl_s1m.canvas.axes.set_title('Se\u00F1al 1 Manipulada')
    ui.mpl_s1m.canvas.axes.plot(t, y1m)
    ui.mpl_s1m.canvas.draw()    
    #genero y grafico la señal 1 manipulada
    if desp_giro:
        y2m = manipulardg(t0, fs, giro, y2)
    else:
        y2m = manipulargd(t0, fs, giro, y2)
    ui.mpl_s2m.canvas.axes.clear()
    ui.mpl_s2m.canvas.axes.set_title('Se\u00F1al 2 Manipulada')
    ui.mpl_s2m.canvas.axes.plot(t, y2m)
    ui.mpl_s2m.canvas.draw()
    
    # generero y grafico las transformadas de fourier de y1
    freqs1 = np.fft.rfftfreq(len(y1), 1/fs)
    fft1 = abs(np.fft.rfft(y1))
    fft1 = fft1 / max(fft1) * amp1
    ui.mpl_fft1.canvas.axes.clear()
    ui.mpl_fft1.canvas.axes.set_title('Se\u00F1al 1 (Frecuencia)')
    ui.mpl_fft1.canvas.axes.plot(freqs1, fft1)
    ui.mpl_fft1.canvas.axes.set_xlim(-1, fs/(per1*25))
    ui.mpl_fft1.canvas.draw()
    # generero y grafico las transformadas de fourier de y2
    freqs2 = np.fft.rfftfreq(len(y2), 1/fs)
    fft2 = abs(np.fft.rfft(y2))
    fft2 = fft2 / max(fft2) * amp2
    ui.mpl_fft2.canvas.axes.clear()
    ui.mpl_fft2.canvas.axes.set_title('Se\u00F1al 2 (Frecuencia)')
    line1, = ui.mpl_fft2.canvas.axes.plot(freqs2, fft2)
    ui.mpl_fft2.canvas.axes.set_xlim(-1, fs/(per2*25))
    ui.mpl_fft2.canvas.draw()
    
    #grafico ambas señales en la vista convolucion
    ui.mpl_both.canvas.axes.clear()
    ui.mpl_both.canvas.axes.set_title('Se\u00F1ales 1 y 2')
    ui.mpl_both.canvas.axes.plot(t, y1, 'r')
    ui.mpl_both.canvas.axes.plot(t, y2, '--g')
    ui.mpl_both.canvas.draw()
    #variables convolucion
        
    conv = signal.convolve(y1, y2)
    conv = conv/max(conv)
    tc = np.linspace(ti*2, ts*2, len(conv))
    zeros = np.zeros(len(tc)-len(t))
    y1 = np.hstack([zeros,y1])
    y2 = np.hstack([y2[::-1],zeros])
    ui.mpl_both.canvas.axes.clear()
    ui.mpl_both.canvas.axes.set_title('Convolucion')
    ui.mpl_both.canvas.axes.plot(tc, conv, 'r')
    ui.mpl_both.canvas.draw()
      
    ui.mpl_convo.canvas.axes.clear()
    ui.mpl_convo.canvas.axes.set_title('Convolucion')
    line1, = ui.mpl_convo.canvas.axes.plot(tc, y2, '--g')
    ui.mpl_convo.canvas.axes.plot(tc, y1, 'b')
    line2, = ui.mpl_convo.canvas.axes.plot(tc, np.zeros(len(tc)), 'r')
    ui.mpl_convo.canvas.draw()
    ventana = int(fs/10)
    for i in range(0, len(tc), ventana):
        desp = np.zeros(i)
        y2b = np.concatenate((desp, y2), axis=None)
        y2b = y2b[:len(y2)]
        conv_update = np.hstack([conv[:i], np.zeros(len(tc))[i:]])
        line1.set_ydata(y2b)
        line2.set_ydata(conv_update)
        ui.mpl_convo.canvas.draw()
        ui.mpl_convo.canvas.flush_events()
        time.sleep(.1)

    
    
    
ui.process.clicked.connect(procesar)



app.exec_()
