
import subprocess
import numpy as np 
import pandas as pd
import scipy.interpolate as scipyint

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# ------------ User Inputs:                                                         ---------------
# ------------ 1. .xml file (compatible w/ Mutation++)                              ---------------
# ------------ 2. table_filename: name of the file the Bprime table will be save to ---------------
# ------------ 3. ts: array of temperatures [K]                                     ---------------
# ------------ 4. ps: array of pressures [Pa]                                       ---------------
# ------------ 5. bs: array of normalized blowing rates                             ---------------
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

xml_file = "tacot26.xml"                    
table_filename = "Bprime_table.bpt"        
ts = np.arange(250,4025,25)                         # Temps from from 250 to 4000 (included) in increments of 25
ps = np.asarray([101.325,1013.25,10132.5,101325])   # Pressures (jumps of a factor of 10) (the jump needs to be equal across all values)
bs = np.arange(-10,10.1,0.1)                        # Blowing/suction rates from -10 to 10 (included) in increments of 0.1

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# ------ This file outputs the Bprime table. The table will contain 5 columns   -------------------
# ------ [Pressure [bar], Bg, Temp [K], Bc, Hw [kJ/Kg]] where Bc is             -------------------
# ------ the ablation rate (related to the surface recession velocity) and Hw   -------------------
# ------ is the formation enthalpy of the chemical equilibrium mixture.         -------------------
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------



# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# ------- No user input required from here on out -------------------------------------------------
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------


bsn = bs[:len(bs)//2+1]                                     # Extract the negative values in bs
bsp = bs[len(bs)//2:]                                       # Extract the positive values in bs


print("Running Mutation++ to generate Table for Bg >= 0 ...")
file = open("bprime_positive.log","w")
proc = subprocess.run(["./generate_bprime_table","-m",xml_file,"-bl","edge","-py","pyro","-T",\
                       "%d:%d:%d"%(ts[0],ts[1]-ts[0],ts[-1]),"-P",\
                       "%1.10f:%1.10f:%1.10f"%(ps[0],ps[1]/ps[0],ps[-1]),"-b",\
                       "%1.5f:%1.5f:%1.5f"%(bsp[0],bsp[1]-bsp[0],bsp[-1])],stdout=file,stderr=subprocess.PIPE,text=True)
file.close()
table_pos = pd.read_csv('bprime_positive.log',header=None,skiprows=1,sep='\s+',engine= 'python').to_numpy()[:,0:5]


print("Running Mutation++ to generate Table for Bg < 0 ...")
file = open("bprime_negative.log","w")
proc = subprocess.run(["./generate_bprime_table","-m",xml_file,"-bl","edge","-py","edge","-T",\
                       "%d:%d:%d"%(ts[0],ts[1]-ts[0],ts[-1]),"-P",\
                       "%1.10f:%1.10f:%1.10f"%(ps[0],ps[1]/ps[0],ps[-1]),"-b",\
                       "%1.5f:%1.5f:%1.5f"%(bsn[0],bsn[1]-bsn[0],bsn[-1])],stdout=file,stderr=subprocess.PIPE,text=True)
file.close()
table_neg = pd.read_csv('bprime_negative.log',header=None,skiprows=1,sep='\s+',engine= 'python').to_numpy()[:,0:5]


ps = ps*1e-5    # Convert pressures to Bar (that's how the tables output the pressure)

# Post process the tables to assemble one unified table

print("Assembling the final B' table ...")
newtable = np.zeros((len(ps)*len(ts)*len(bs),5))

for pk in range (len(ps)):
    idcespk = np.arange(pk*len(bsn)*len(ts),(pk+1)*len(ts)*len(bsn),1)
    tab_pos = table_pos[idcespk,]
    tab_neg = table_neg[idcespk,]
    
    for tk in range (len(ts)):
    
        idcestk = np.argwhere(np.abs(tab_pos[:,2] - ts[tk]) <= 1e-10)[:,-1]
        tab = np.concatenate((tab_neg[idcestk[:-1],],tab_pos[idcestk,]),axis=0)
        
        # Positive values
        bg_pos = tab[len(bs)//2:,1] 
        bf_pos = tab[len(bs)//2:,1] + tab[len(bs)//2:,-2]
        
        # Negative values
        bf_neg = tab[0:len(bs)//2,1]
        bg_neg = tab[0:len(bs)//2,1] - tab[0:len(bs)//2,-2]
        
        bg = np.concatenate((bg_neg.reshape(-1,1),bg_pos.reshape(-1,1))).reshape(-1)
        bf = np.concatenate((bf_neg.reshape(-1,1),bf_pos.reshape(-1,1))).reshape(-1)
        
        bc = tab[:,-2]
        hw = tab[:,-1]
        
        fbc = scipyint.interp1d(bg,bc,kind='linear',fill_value='extrapolate')
        fhw = scipyint.interp1d(bg,hw,kind='linear',fill_value='extrapolate')
        
        bg = bs.copy()
        bc = fbc(bg)
        hw = fhw(bg)
        
        idces = pk*len(bs)*len(ts) + np.arange(0,len(ts)*len(bs),len(ts)) + tk
        
        newtable[idces,0] = ps[pk]
        newtable[idces,1] = bg
        newtable[idces,2] = ts[tk]
        newtable[idces,3] = bc
        newtable[idces,4] = hw*1e3
        

# Save the table to file
print("Saving the Bprime table to file %s"%(table_filename))
file = open(table_filename,"w")
file.write("// P (bar) Bg T (K) Bc Hw(kJ/kg)\n")
for k in range (len(ps)):
    for j in range (len(bs)):    
        for i in range (len(ts)):
            idx = i + j*len(ts) + k*len(bs)*len(ts)
            vec = newtable[idx,]
            for zz in range (len(vec)):
                file.write("%1.10E  "%vec[zz])
            file.write("\n")
        file.write("\n")
file.close()


# Wipe outputs of Mutation++
proc = subprocess.run(["rm","bprime_positive.log"])
proc = subprocess.run(["rm","bprime_negative.log"])  
        






























