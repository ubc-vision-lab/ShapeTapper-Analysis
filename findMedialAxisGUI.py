#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 13:55:18 2018

@author: visionlab
"""

from Tkinter import *
from tkMessageBox import *
import tkFileDialog
import os
import mat_analysis as ma


def callback_dir():
    
    showinfo('Data', 'Please select directory containing shape analysis files (.png & .mat)')
    img_path = tkFileDialog.askdirectory()
    img_names = [fname for fname in os.listdir(img_path) if fname.endswith(".png")]
    img_names.sort()
    img_mats = [fname for fname in os.listdir(img_path) if fname.endswith(".mat")]
    img_mats.sort()
    flist = [f[:len(f)-7] for f in img_names]
    
    showinfo('Data', 'Please select directory containing touch data files (.mat)')
    data_path = tkFileDialog.askdirectory()
    data_mats = [fname for fname in os.listdir(data_path) if fname.endswith(".mat")]
    data_mats.sort()
    
    if len(img_names) == len(data_mats) and len(img_names) == len(img_mats):
        confirm_msg = "You have selected the following files:\n\n {}\n\nSelect output directory & continue.".format(flist)
        if askyesno('Verify', confirm_msg):
            out_path = tkFileDialog.askdirectory()
            for i,img in enumerate(img_names):
                ma.matAnalysis(img_names[i], img_mats[i], img_path, data_mats[i], data_path, out_path)
            showwarning('Success', 'Done!')
    else:
        showinfo('Error', 'Number of files do not match!')


root = Tk()
root.title("Medial Axis Analysis")
w = root.winfo_reqwidth()
h = root.winfo_reqheight()
ws = root.winfo_screenwidth()
hs = root.winfo_screenheight()
x = (ws/2) - (w/2)
y = (hs/2) - (h/2)
root.geometry('+%d+%d' % (x, y))
root.geometry("200x200")
errmsg = 'Error!'
Button(text='Directory Open', command=callback_dir).grid(row=0, sticky=N+S+E+W)
Button(text='Quit', command=root.destroy).grid(row=1, sticky=N+S+E+W)
mainloop()