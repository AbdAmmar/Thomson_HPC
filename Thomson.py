import customtkinter
import tkinter
from numpy import pi
from math import sqrt
from numpy import random
import numpy as np
from tkinter.filedialog import asksaveasfile
import os , sys


customtkinter.set_appearance_mode("System")  # Modes: "System" (standard), "Dark", "Light"
customtkinter.set_default_color_theme("green")  # Themes: "blue" (standard), "green", "dark-blue"


class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        self.geometry("1020x760")
        self.title("Thomson Input generator")
        self.minsize(1020, 760)

        #self.grid_rowconfigure(0, weight=1)

        self.textbox = customtkinter.CTkTextbox(master=self , width = 500)
        self.textbox.grid(row=0, column=1, padx=20, pady=(40, 20) , sticky="ns")
        self.textbox.configure(state='disabled')
        
        self.tabview=customtkinter.CTkTabview(master=self)
        self.tabview.grid(row=0,column=0,padx=20, columnspan=1 , pady=(20,20), sticky="ns")
        
        left = self.tabview.add("Custom")  
        right = self.tabview.add("Preload")  
        self.tabview.set("Custom")
      
        
        self.continer1 = customtkinter.CTkFrame(master=left)
        self.continer1.grid(row = 0 , column=0, padx=20, columnspan=1 , pady=(20,20), sticky="ns")
        
        self.continer3 = customtkinter.CTkFrame(master=right)
        self.continer3.grid(row = 0 , column=0, padx=20, columnspan=1 , pady=(20,20), sticky="ns")
        
        ####################################  tab one  ######################################################
        
        self.label_Di = customtkinter.CTkLabel(master=self.continer1 , text="Dimension :"  )
        self.label_Di.grid(row=0 ,column=0, padx=10, pady=5 ,sticky="nsew" )
        
        self.combobox = customtkinter.CTkComboBox(master=self.continer1, values=["1", "2" , "3"])
        self.combobox.grid(row=0, column=1, padx=50, pady=5 , sticky="w" )
        self.combobox.configure(width = 100 , justify = "center" )
        self.combobox.set("3")
        
        self.label_Di = customtkinter.CTkLabel(master=self.continer1 , text="  Number  of electrons :"  )
        self.label_Di.grid(row=1 ,column=0, padx=10, pady=5, sticky="nsew" )
        
        self.electron_N = customtkinter.CTkTextbox(master=self.continer1)
        self.electron_N.grid(row=1 , column=1 , padx=10, pady=5)
        self.electron_N.configure(width = 60, height = 20 )
        self.electron_N.insert("0.0", "10")
        
        self.label_tolerance = customtkinter.CTkLabel(master=self.continer1 , text=" Tolerance :"  )
        self.label_tolerance.grid(row=2 ,column=0, padx=10, pady=5, sticky="nsew" )
        
        
        self.tolerance_N = customtkinter.CTkTextbox(master=self.continer1)
        self.tolerance_N.grid(row=2 , column=1 , padx=10, pady=5)
        self.tolerance_N.configure(width = 60, height = 20 )
        self.tolerance_N.insert("0.0", "1e-8")
        
        
        self.label_itermax = customtkinter.CTkLabel(master=self.continer1 , text=" Total number of iteration :"  )
        self.label_itermax.grid(row=3 ,column=0, padx=10, pady=5, sticky="nsew" )
        
        
        self.itermax_N = customtkinter.CTkTextbox(master=self.continer1)
        self.itermax_N.grid(row=3 , column=1 , padx=10, pady=5)
        self.itermax_N.configure(width = 60, height = 20 )
        self.itermax_N.insert("0.0", "10000")
        
        
        ############################### the size of the box #################################
        
        self.continer2 = customtkinter.CTkFrame(master=self.continer1)
        self.continer2.grid(row = 4 , column=0 ,columnspan=2 ,  padx=(20,20), pady=(20,20), sticky="we")
    
        
        self.label_size = customtkinter.CTkLabel(master=self.continer2 , text=" The Box Size"  )
        self.label_size.grid(row=0 ,column=0,  padx=(10,10), pady=5 )
        
        
        self.label_Lx = customtkinter.CTkLabel(master=self.continer2 , text=" Lx"  )
        self.label_Lx.grid(row=1 ,column=0,  padx=(10,10), pady=5 )
        
        self.Lx = customtkinter.CTkTextbox(master=self.continer2)
        self.Lx.grid(row=1 ,column=1, padx=10, pady=5)
        self.Lx.configure(width = 150, height = 20 )
        self.Lx.insert("0.0", 2*pi)
        
        self.label_Ly = customtkinter.CTkLabel(master=self.continer2 , text=" Ly"  )
        self.label_Ly.grid(row=2 ,column=0,  padx=(10,10), pady=5 )
        
        self.Ly = customtkinter.CTkTextbox(master=self.continer2)
        self.Ly.grid(row=2 ,column=1, padx=10, pady=5)
        self.Ly.configure(width = 150, height = 20 )
        self.Ly.insert("0.0", 2*pi)
        
        
        self.label_Lz = customtkinter.CTkLabel(master=self.continer2 , text=" Lz"  )
        self.label_Lz.grid(row=3 ,column=0,  padx=(10,10), pady=5 )
        
        self.Lz = customtkinter.CTkTextbox(master=self.continer2)
        self.Lz.grid(row=3 ,column=1, padx=10, pady=(5,20))
        self.Lz.configure(width = 150, height = 20 )        
        self.Lz.insert("0.0", 2*pi)
       
       ##########################################################################################
       
       
       
       ###################################### check box #########################################
        
        self.checkbox_R = customtkinter.CTkCheckBox(master=self.continer1, text="Random Geometry", onvalue="random", offvalue="")
        self.checkbox_R.grid(row=5 ,column=1, padx=10, pady=(0,10) , sticky="w"  )
        
        self.checkbox_A = customtkinter.CTkCheckBox(master=self.continer1, text="Animation", onvalue="animation", offvalue="")
        self.checkbox_A.grid(row=5 ,column=0, padx=10, pady=(0,10) , sticky="w"  )
        
        self.checkbox_H = customtkinter.CTkCheckBox(master=self.continer1, text="Show Hessian Matrix", onvalue="hessian", offvalue="")
        self.checkbox_H.grid(row=6 ,column=0, padx=10, pady=10 , sticky="w"  )
        
        self.checkbox_S = customtkinter.CTkCheckBox(master=self.continer1, text="Show All Result", onvalue="show", offvalue="")
        self.checkbox_S.grid(row=7 ,column=0, padx=10, pady=10 , sticky="w"  )
        
        self.checkbox_M = customtkinter.CTkCheckBox(master=self.continer1, text="Multiply by the Box size", onvalue="multiply", offvalue="")
        self.checkbox_M.grid(row=6 ,column=1,columnspan=3,  padx=10, pady=10 , sticky="w"  )
        
        
        ############################################################################################
        
        
        self.button_write = customtkinter.CTkButton(master=self.continer1, command=self.button_write_callback, text="Write")
        self.button_write.grid(row=9, column=0 ,  padx=20, pady=(20,20), sticky="nsew")
        
        self.button_clear = customtkinter.CTkButton(master=self.continer1, command=self.button_clear_callback, text="Clear")
        self.button_clear.grid(row=9, column=1 ,  padx=20, pady=(20,20), sticky="nsew")
        
        self.button_save = customtkinter.CTkButton(master=self.continer1, command=self.button_save_callback, text="Save")
        self.button_save.grid(row=11, column=0, columnspan=3 ,  padx=20, pady=(0,20), sticky="nsew")
        
        self.button_run = customtkinter.CTkButton(master=self.continer1, command=self.button_run_callback, text="Run")
        self.button_run.grid(row=12, column=0, columnspan=3 ,  padx=20, pady=(0,20), sticky="nsew")
        
        
        self.tabview=customtkinter.CTkTabview(master=self)
        
        
       ################################################ tab two ##################################
       
        self.label_CS = customtkinter.CTkLabel(master=self.continer3 , text="Crystal Structure "  )
        self.label_CS.grid(row=0 ,column=0, padx=10, pady=5 ,sticky="nsew" )
        
        self.combobox_CS = customtkinter.CTkComboBox(master=self.continer3, values=["FCC (3D)", "BCC (3D)" , "SC (3D)" ,"HEX (2D)" , "SL (2D)" , "ED (1D)"])
        self.combobox_CS.grid(row=0, column=1, padx=50, pady=5 , sticky="w" )
        self.combobox_CS.configure(width = 100 , justify = "center" )
        self.combobox_CS.set("FCC (3D)")
        
        
        self.label_CS_EN = customtkinter.CTkLabel(master=self.continer3 , text="Unit-cell per dimension"  )
        self.label_CS_EN.grid(row=1 ,column=0, padx=10, pady=5 ,sticky="nsew" )
        
        self.electron_N_EN = customtkinter.CTkTextbox(master=self.continer3)
        self.electron_N_EN.grid(row=1 , column=1 , padx=10, pady=5)
        self.electron_N_EN.configure(width = 60, height = 20 )
        self.electron_N_EN.insert("0.0", "1")
        
        self.label_tolerance_EN = customtkinter.CTkLabel(master=self.continer3 , text=" Tolerance :"  )
        self.label_tolerance_EN.grid(row=2 ,column=0, padx=10, pady=5, sticky="nsew" )
        
        
        self.tolerance_N_EN = customtkinter.CTkTextbox(master=self.continer3)
        self.tolerance_N_EN.grid(row=2 , column=1 , padx=10, pady=5)
        self.tolerance_N_EN.configure(width = 60, height = 20 )
        self.tolerance_N_EN.insert("0.0", "1e-8")
        
        
        self.label_itermax_EN = customtkinter.CTkLabel(master=self.continer3 , text=" Total number of iteration :"  )
        self.label_itermax_EN.grid(row=3 ,column=0, padx=10, pady=5, sticky="nsew" )
        
        
        self.itermax_N_EN = customtkinter.CTkTextbox(master=self.continer3)
        self.itermax_N_EN.grid(row=3 , column=1 , padx=10, pady=5)
        self.itermax_N_EN.configure(width = 60, height = 20 )
        self.itermax_N_EN.insert("0.0", "10000")
        
        self.continer4 = customtkinter.CTkFrame(master=self.continer3)
        self.continer4.grid(row = 4 , column=0 ,columnspan=2 ,  padx=(20,20), pady=(20,20), sticky="we")
    
        
        self.label_size = customtkinter.CTkLabel(master=self.continer4 , text=" The Box Size"  )
        self.label_size.grid(row=0 ,column=0,  padx=(10,10), pady=5 )
        
        
        self.label_Lx_EN = customtkinter.CTkLabel(master=self.continer4 , text=" Lx"  )
        self.label_Lx_EN.grid(row=1 ,column=0,  padx=(10,10), pady=5 )
        
        self.Lx_EN = customtkinter.CTkTextbox(master=self.continer4)
        self.Lx_EN.grid(row=1 ,column=1, padx=10, pady=5)
        self.Lx_EN.configure(width = 150, height = 20 )
        self.Lx_EN.insert("0.0", 2*pi)
        
        self.label_Ly_EN = customtkinter.CTkLabel(master=self.continer4 , text=" Ly"  )
        self.label_Ly_EN.grid(row=2 ,column=0,  padx=(10,10), pady=5 )
        
        self.Ly_EN = customtkinter.CTkTextbox(master=self.continer4)
        self.Ly_EN.grid(row=2 ,column=1, padx=10, pady=5)
        self.Ly_EN.configure(width = 150, height = 20 )
        self.Ly_EN.insert("0.0", 2*pi)
        
        
        self.label_Lz_EN = customtkinter.CTkLabel(master=self.continer4 , text=" Lz"  )
        self.label_Lz_EN.grid(row=3 ,column=0,  padx=(10,10), pady=5 )
        
        self.Lz_EN = customtkinter.CTkTextbox(master=self.continer4)
        self.Lz_EN.grid(row=3 ,column=1, padx=10, pady=(5,20))
        self.Lz_EN.configure(width = 150, height = 20 )        
        self.Lz_EN.insert("0.0", 2*pi)
        
        self.checkbox_A_EN = customtkinter.CTkCheckBox(master=self.continer3, text="Animation", onvalue="animation", offvalue="")
        self.checkbox_A_EN.grid(row=5 ,column=0, padx=10, pady=(0,10) , sticky="w"  )
        
        self.checkbox_H_EN = customtkinter.CTkCheckBox(master=self.continer3, text="Show Hessian Matrix", onvalue="hessian", offvalue="")
        self.checkbox_H_EN.grid(row=6 ,column=0, padx=10, pady=10 , sticky="w"  )
        
        self.checkbox_S_EN = customtkinter.CTkCheckBox(master=self.continer3, text="Show All Result", onvalue="show", offvalue="")
        self.checkbox_S_EN.grid(row=7 ,column=0, padx=10, pady=10 , sticky="w"  )
        
        self.button_write_EN = customtkinter.CTkButton(master=self.continer3, command=self.button_write_EN_callback, text="Write")
        self.button_write_EN.grid(row=9, column=0 ,  padx=20, pady=(20,20), sticky="nsew")
        
        self.button_clear_EN = customtkinter.CTkButton(master=self.continer3, command=self.button_clear_EN_callback, text="Clear")
        self.button_clear_EN.grid(row=9, column=1 ,  padx=20, pady=(20,20), sticky="nsew")
        
        self.button_save_EN = customtkinter.CTkButton(master=self.continer3, command=self.button_save_EN_callback, text="Save")
        self.button_save_EN.grid(row=10, column=0, columnspan=3 ,  padx=20, pady=(0,20), sticky="nsew")
        
        self.button_run_EN = customtkinter.CTkButton(master=self.continer3, command=self.button_run_EN_callback, text="Run")
        self.button_run_EN.grid(row=11, column=0, columnspan=3 ,  padx=20, pady=(0,20), sticky="nsew")
        


        self.Lx_EN.configure(state='disabled')
        self.Ly_EN.configure(state='disabled')
        self.Lz_EN.configure(state='disabled')

        ###########################################################################################
        
        
        ###################### functions ########################
        

    def button_write_callback(self):
        self.textbox.configure(state='normal')
        self.textbox.insert("insert","dimension:  " + self.combobox.get() + "\n")
        self.textbox.insert("insert","electron:  " + self.electron_N.get("0.0", "end") )
        self.textbox.insert("insert","tolerance:  " + self.tolerance_N.get("0.0", "end") )
        self.textbox.insert("insert","itermax:  " + self.itermax_N.get("0.0", "end"))
        self.textbox.insert("insert","box:  " + self.Lx.get("0.0", "end-1c")  +" "+ self.Ly.get("0.0", "end-1c") +"  "+ self.Lz.get("0.0", "end"))
        self.textbox.insert("insert","periodic" + "\n")
        
        
        if self.checkbox_R.get() == "random":
            self.textbox.insert("insert", self.checkbox_R.get() + "\n")   # type: ignore
        if self.checkbox_M.get() == "multiply":
            self.textbox.insert("insert", self.checkbox_M.get() + "\n")   # type: ignore
        if self.checkbox_S.get() == "show":
            self.textbox.insert("insert", self.checkbox_S.get() + "\n")   # type: ignore  
        if self.checkbox_A.get() == "animation":
            self.textbox.insert("insert", self.checkbox_A.get() + "\n")   # type: ignore
        if self.checkbox_H.get() == "hessian":
            self.textbox.insert("insert", self.checkbox_H.get() + "\n")   # type: ignore
        
        if self.checkbox_R.get() == "random":
            self.textbox.insert("insert","geometry" + "\n")
            x = random.rand(int(self.electron_N.get("0.0", "end-1c")))
            y = random.rand(int(self.electron_N.get("0.0", "end-1c")))
            z = random.rand(int(self.electron_N.get("0.0", "end-1c")))
            if str(self.combobox.get()) == "3":
                for i in range(int(self.electron_N.get("0.0", "end-1c"))):
                    self.textbox.insert("insert", "  " + str(i+1) + "     " + str(x[i]).replace('[','  ').replace(']','  ') + "    " + str(y[i]).replace('[','  ').replace(']','  ') + "    " + str(z[i]).replace('[','  ').replace(']','  ') + "\n")  # type: ignore
            elif str(self.combobox.get()) == "2":
                for i in range(int(self.electron_N.get("0.0", "end-1c"))):
                    self.textbox.insert("insert", "  " + str(i+1) + "     " + str(x[i]).replace('[','  ').replace(']','  ') + "    " + str(y[i]).replace('[','  ').replace(']','  ') + "\n")  # type: ignore                
            elif str(self.combobox.get()) == "1":
                for i in range(int(self.electron_N.get("0.0", "end-1c"))):
                    self.textbox.insert("insert", "  " + str(i+1) + "     " + str(x[i]).replace('[','  ').replace(']','  ') +  "\n")  # type: ignore                
        else:
            self.textbox.insert("insert","# Change the geometry as you want" + "\n")
            self.textbox.insert("insert","geometry" + "\n")
            x = np.zeros(int(self.electron_N.get("0.0", "end-1c")))
            y = np.zeros(int(self.electron_N.get("0.0", "end-1c")))
            z = np.zeros(int(self.electron_N.get("0.0", "end-1c")))
            if str(self.combobox.get()) == "3":        
                for i in range(int(self.electron_N.get("0.0", "end-1c"))):
                    self.textbox.insert("insert", "  " + str(i+1) + "    " + str(x[i]).replace('[','  ').replace(']','  ') + "    " + str(y[i]).replace('[','  ').replace(']','  ') + "    " + str(z[i]).replace('[','  ').replace(']','  ') + "\n")
            elif str(self.combobox.get()) == "2":        
                for i in range(int(self.electron_N.get("0.0", "end-1c"))):
                    self.textbox.insert("insert", "  " + str(i+1) + "    " + str(x[i]).replace('[','  ').replace(']','  ') + "    " + str(y[i]).replace('[','  ').replace(']','  ') +  "\n")
            elif str(self.combobox.get()) == "1":        
                for i in range(int(self.electron_N.get("0.0", "end-1c"))):
                    self.textbox.insert("insert", "  " + str(i+1) + "    " + str(x[i]).replace('[','  ').replace(']','  ') +  "\n")
        
        self.textbox.insert("insert", "\n")
        self.textbox.insert("insert", "______________________________________________________" + "\n")
        self.textbox.insert("insert", "\n")
        self.textbox.configure(state='disabled')



    def button_write_EN_callback(self):
        self.textbox.configure(state='normal')
        
        if self.combobox_CS.get() == "FCC (3D)" or self.combobox_CS.get() == "BCC (3D)" or self.combobox_CS.get() == "SC (3D)":
            self.textbox.insert("insert","dimension:  3 "  + "\n")
        elif self.combobox_CS.get() == "SL (2D)" or self.combobox_CS.get() == "HEX (2D)":
            self.textbox.insert("insert","dimension:  2 "  + "\n")
        elif self.combobox_CS.get() == "ED (1D)":
            self.textbox.insert("insert","dimension:  1 "  + "\n")
        
        if self.combobox_CS.get() == "FCC (3D)":
            self.Lx_EN.configure(state='normal')
            self.Ly_EN.configure(state='normal')
            self.Lz_EN.configure(state='normal')
            self.Ly_EN.delete("0.0","end")
            self.Ly_EN.insert("0.0", 2*pi)
            self.Lx_EN.configure(state='disabled')
            self.Ly_EN.configure(state='disabled')
            self.Lz_EN.configure(state='disabled')
            x=4*int(self.electron_N_EN.get("0.0", "end"))**3
        elif self.combobox_CS.get() == "BCC (3D)":
            self.Lx_EN.configure(state='normal')
            self.Ly_EN.configure(state='normal')
            self.Lz_EN.configure(state='normal')
            self.Ly_EN.delete("0.0","end")
            self.Ly_EN.insert("0.0", 2*pi)
            self.Lx_EN.configure(state='disabled')
            self.Ly_EN.configure(state='disabled')
            self.Lz_EN.configure(state='disabled')
            x=2*int(self.electron_N_EN.get("0.0", "end"))**3
        elif self.combobox_CS.get() == "SC (3D)":
            self.Lx_EN.configure(state='normal')
            self.Ly_EN.configure(state='normal')
            self.Lz_EN.configure(state='normal')
            self.Ly_EN.delete("0.0","end")
            self.Ly_EN.insert("0.0", 2*pi)
            self.Lx_EN.configure(state='disabled')
            self.Ly_EN.configure(state='disabled')
            self.Lz_EN.configure(state='disabled')
            x=int(self.electron_N_EN.get("0.0", "end"))**3
        elif self.combobox_CS.get() == "HEX (2D)":
            self.Lx_EN.configure(state='normal')
            self.Ly_EN.configure(state='normal')
            self.Lz_EN.configure(state='normal')
            self.Ly_EN.delete("0.0","end")
            self.Ly_EN.insert("0.0", sqrt(3)*pi)
            self.Lx_EN.configure(state='disabled')
            self.Ly_EN.configure(state='disabled')
            self.Lz_EN.configure(state='disabled')
            x=4*int(self.electron_N_EN.get("0.0", "end"))**2
        elif self.combobox_CS.get() == "SL (2D)":
            self.Lx_EN.configure(state='normal')
            self.Ly_EN.configure(state='normal')
            self.Lz_EN.configure(state='normal')
            self.Ly_EN.delete("0.0","end")
            self.Ly_EN.insert("0.0", 2*pi)
            self.Lx_EN.configure(state='disabled')
            self.Ly_EN.configure(state='disabled')
            self.Lz_EN.configure(state='disabled')
            x=int(self.electron_N_EN.get("0.0", "end"))**2
        elif self.combobox_CS.get() == "ED (1D)":
            self.Lx_EN.configure(state='normal')
            self.Ly_EN.configure(state='normal')
            self.Lz_EN.configure(state='normal')
            self.Ly_EN.delete("0.0","end")
            self.Ly_EN.insert("0.0", 2*pi)
            self.Lx_EN.configure(state='disabled')
            self.Ly_EN.configure(state='disabled')
            self.Lz_EN.configure(state='disabled')
            x=int(self.electron_N_EN.get("0.0", "end"))
            
        self.textbox.insert("insert","electron:  " + str(x)+"\n") # type: ignore
        self.textbox.insert("insert","tolerance:  " + self.tolerance_N_EN.get("0.0", "end") )
        self.textbox.insert("insert","itermax:  " + self.itermax_N_EN.get("0.0", "end"))
        self.textbox.insert("insert","box:  " + self.Lx_EN.get("0.0", "end-1c")  +"  "+ self.Ly_EN.get("0.0", "end-1c") +"  "+ self.Lz_EN.get("0.0", "end"))
        self.textbox.insert("insert","periodic" + "\n")
        
        if self.checkbox_S_EN.get() == "show":
            self.textbox.insert("insert", self.checkbox_S_EN.get() + "\n")   # type: ignore  
        if self.checkbox_A_EN.get() == "animation":
            self.textbox.insert("insert", self.checkbox_A_EN.get() + "\n")   # type: ignore
        if self.checkbox_H_EN.get() == "hessian":
            self.textbox.insert("insert", self.checkbox_H_EN.get() + "\n")   # type: ignore
        
           
        
        if self.combobox_CS.get() == "FCC (3D)":
            d1 = 1/int(self.electron_N_EN.get("0.0", "end"))
            d2 = d1/2
            icount = 0
            
            self.textbox.insert("insert", "multiply" + "\n")   # type: ignore
            self.textbox.insert("insert","geometry"  + "\n") 
            
            
            for i in range(int(self.electron_N_EN.get("0.0", "end"))):
                for j in range(int(self.electron_N_EN.get("0.0", "end"))):
                    for k in range(int(self.electron_N_EN.get("0.0", "end"))):
                        icount += 1
                        self.textbox.insert("insert", str(icount) + "\t" + str(round(d1*i,18)) + "\t \t" + str(round(d1*j,18)) + "\t \t" + str(round(d1*k,18)) + "\n")
                        
            for i in range(int(self.electron_N_EN.get("0.0", "end"))):
                for j in range(int(self.electron_N_EN.get("0.0", "end"))):
                    for k in range(int(self.electron_N_EN.get("0.0", "end"))):
                        icount += 1
                        self.textbox.insert("insert", str(icount) + "\t" + str(round(d1*i+d2,18)) + "\t \t" + str(round(d1*j+d2,18)) + "\t \t" + str(round(d1*k,18)) + "\n")
            
            for i in range(int(self.electron_N_EN.get("0.0", "end"))):
                for j in range(int(self.electron_N_EN.get("0.0", "end"))):
                    for k in range(int(self.electron_N_EN.get("0.0", "end"))):
                        icount += 1
                        self.textbox.insert("insert", str(icount) + "\t" + str(round(d1*i,18)) + "\t \t" + str(round(d1*j+d2,18)) + "\t \t" + str(round(d1*k+d2,18)) + "\n")
                        
            for i in range(int(self.electron_N_EN.get("0.0", "end"))):
                for j in range(int(self.electron_N_EN.get("0.0", "end"))):
                    for k in range(int(self.electron_N_EN.get("0.0", "end"))):
                        icount += 1
                        self.textbox.insert("insert", str(icount) + "\t" + str(round(d1*i+d2,18)) + "\t \t" + str(round(d1*j,18)) + "\t \t" + str(round(d1*k+d2,18)) + "\n")
                        
                        
                        
        if self.combobox_CS.get() == "BCC (3D)":
            d1 = 1/int(self.electron_N_EN.get("0.0", "end"))
            d2 = d1/2
            icount = 0
            
            self.textbox.insert("insert", "multiply" + "\n")   # type: ignore
            self.textbox.insert("insert","geometry" + "\n") 
            
            for i in range(int(self.electron_N_EN.get("0.0", "end"))):
                for j in range(int(self.electron_N_EN.get("0.0", "end"))):
                    for k in range(int(self.electron_N_EN.get("0.0", "end"))):
                        icount += 1
                        self.textbox.insert("insert", str(icount) + "\t" + str(round(d1*i,18)) + "\t \t" + str(round(d1*j,18)) + "\t \t" + str(round(d1*k,18)) + "\n") 
            
            for i in range(int(self.electron_N_EN.get("0.0", "end"))):
                for j in range(int(self.electron_N_EN.get("0.0", "end"))):
                    for k in range(int(self.electron_N_EN.get("0.0", "end"))):
                        icount += 1
                        self.textbox.insert("insert", str(icount) + "\t" + str(round(d1*i+d2,18)) + "\t \t" + str(round(d1*j+d2,18)) + "\t \t" + str(round(d1*k+d2,18)) + "\n")               
        
            
        if self.combobox_CS.get() == "SC (3D)":
            d1 = 1/int(self.electron_N_EN.get("0.0", "end"))
            icount = 0
            
            self.textbox.insert("insert", "multiply" + "\n")   # type: ignore
            self.textbox.insert("insert","geometry" + "\n") 
            
            for i in range(int(self.electron_N_EN.get("0.0", "end"))):
                for j in range(int(self.electron_N_EN.get("0.0", "end"))):
                    for k in range(int(self.electron_N_EN.get("0.0", "end"))):
                        icount += 1
                        self.textbox.insert("insert", str(icount) + "\t" + str(round(d1*i,18)) + "\t \t" + str(round(d1*j,18)) + "\t \t" + str(round(d1*k,18)) + "\n") 
        
        if self.combobox_CS.get() == "SL (2D)":
            d1 = 1/int(self.electron_N_EN.get("0.0", "end"))
            icount = 0
            
            self.textbox.insert("insert", "multiply" + "\n")   # type: ignore
            self.textbox.insert("insert","geometry" + "\n") 
            
            for i in range(int(self.electron_N_EN.get("0.0", "end"))):
                for j in range(int(self.electron_N_EN.get("0.0", "end"))):
                    icount += 1
                    self.textbox.insert("insert", str(icount) + "\t" + str(round(d1*i,20)) + "\t \t \t" + str(round(d1*j,20)) + "\t" + "\n") 
                    
        if self.combobox_CS.get() == "HEX (2D)":
            dx = pi/int(self.electron_N_EN.get("0.0", "end"))
            dy = dx*sqrt(3)/2
            icount = 0

            self.textbox.insert("insert", "multiply" + "\n")   # type: ignore

            self.textbox.insert("insert","geometry" + "\n") 
            for j in range(int(self.electron_N_EN.get("0.0", "end"))):
                
                for i in range(2*int(self.electron_N_EN.get("0.0", "end"))):
                    icount += 1
                    self.textbox.insert("insert", str(icount) + "\t" + str(round(dx*i/(2*pi),20)) + "\t \t \t" + str(round(2*dy*j/(sqrt(3)*pi),20)) + "\t" + "\n") 
                
                for i in range(2*int(self.electron_N_EN.get("0.0", "end"))):
                    icount += 1
                    self.textbox.insert("insert", str(icount) + "\t" + str(round((dx*(i+0.5)/(2*pi)),20)) + "\t \t \t" + str(round(((2*j+1)*dy/(sqrt(3)*pi)),20)) + "\t" + "\n")
                                       
        
        if self.combobox_CS.get() == "ED (1D)":
            dx = 1/int(self.electron_N_EN.get("0.0", "end"))
            icount = 0
            self.textbox.insert("insert", "multiply" + "\n")   # type: ignore
            self.textbox.insert("insert","geometry" + "\n") 
            for i in range(1,int(self.electron_N_EN.get("0.0", "end"))+1):
                icount += 1
                self.textbox.insert("insert", str(icount) + "\t" + str(round(dx*i,20))+ "\n") 
        
        
        
        self.textbox.insert("insert", "\n")
        self.textbox.insert("insert", "______________________________________________________" + "\n")
        self.textbox.insert("insert", "\n")
        self.textbox.configure(state='disabled')
        
        



    def button_clear_callback(self):
        self.textbox.configure(state='normal')
        self.textbox.delete("0.0","end")
        self.textbox.configure(state='disabled')
        
    def button_clear_EN_callback(self):
        self.textbox.configure(state='normal')
        self.textbox.delete("0.0","end")
        self.textbox.configure(state='disabled')

    def button_save_callback(self):
        self.dialog = customtkinter.CTkInputDialog(text="Type The name without extension:", title="Input file save")
        
        file = open(str(self.dialog.get_input()) + ".inp", "w")  
        file.write(self.textbox.get("0.0","end-1c"))  
        file.close()
    
    def button_save_EN_callback(self):
        self.dialog = customtkinter.CTkInputDialog(text="Type The name without extension:", title="Input file save")
        
        file = open(str(self.dialog.get_input()) + ".inp", "w")  
        file.write(self.textbox.get("0.0","end-1c"))  
        file.close()    

    def button_run_callback(self):
        file = open("null.inp","w")
        file.write(self.textbox.get("0.0","end-1c"))
        file.close()
        os.system("./Thomson null.inp")
        os.system("rm -rf null.inp")
        os.system("echo '\n'")
        os.system("echo    '((Finished, you could continue for another calculation ))'   ")

        
        
    def button_run_EN_callback(self):
        file = open("null.inp","w")
        file.write(self.textbox.get("0.0","end-1c"))
        file.close()
        os.system("./Thomson null.inp")
        os.system("rm -rf null.inp")
        os.system("echo '\n'")
        os.system("echo    '((Finished, you could continue for another calculation ))'   ")
    

if __name__ == "__main__":
    app = App()
    app.mainloop()