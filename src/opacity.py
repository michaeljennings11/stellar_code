import numpy as np
from scipy.interpolate import RegularGridInterpolator, NearestNDInterpolator, LinearNDInterpolator
from io import StringIO

opal_file = './opacities/GN93/GN93hz'

tables_begin = 240
tables_head = 0
table_end = 77
table_begin = 0
table_length = table_end - table_begin

def open_opal(file_path):
    with open(file_path,'r') as f:
        flines = f.readlines()
    tables_begin = 240
    tables_head = 0
    table_lines = len(flines) - tables_begin - tables_head
    table_end = 77
    table_begin = 0
    table_length = table_end - table_begin
    table_num   = int(table_lines / table_length)
    
    summary_begin = 62
    summary_end = 188
    summary_text = "".join(flines[summary_begin:summary_end])
    
    
    sconv = lambda s : float(s.split("=")[-1]) 
    tconv = lambda s : float(s.split("#")[-1])   
    myConverters = {}
    myConverters[0] = tconv
    myConverters[2] = sconv
    myConverters[3] = sconv
    myConverters[4] = sconv
    myConverters[5] = sconv
    myConverters[6] = sconv
    myConverters
    
    summary = np.genfromtxt(StringIO(summary_text),
        usecols=[0,2,3,4,5,6],
        encoding=None,
        comments="*",
        delimiter=[11,16,12,9,9,11,11],
        autostrip=True,
        converters=myConverters
    )
    
    comp = np.empty((table_num,5))
    
    for tn in range(126):
        X = summary[tn][1]
        Y = summary[tn][2]
        Z = summary[tn][3]
        dXc = summary[tn][4]
        dXo = summary[tn][5]
        assert X + Y + Z - 1.0 < 1e-10, "X: %g, Y: %g, Z: %g S: %g" % (X,Y,Z,X+Y+Z)

        comp[tn] = np.array([X,Y,dXc,dXo,tn])
        
    return comp,flines

def find_table(point,comps):
    points = comps[:,:4]
    values = comps[:,4]
    tables_number = NearestNDInterpolator(points,values)
    return int(tables_number(point))+1

def interp_opac(point,file_path,method="regular",return_values=False):
    comps,flines = open_opal(file_path)
    tab_num = find_table(point,comps)
    
    tables_begin = 240
    tables_head = 0
    table_lines = len(flines) - tables_begin - tables_head
    table_end = 77
    table_begin = 0
    table_length = table_end - table_begin
    table_num   = table_lines / table_length
    assert table_num == 126
    table_num = int(table_num)
    
    kappas = np.empty((70,19))
    kappas[:] = np.nan
    sheader = np.array(flines[tables_begin+(table_length*(tab_num-1))+5:tables_begin+(table_length*(tab_num-1))+5+1][0].split()[1:])
    logRs = [float(s) for s in sheader]
    rows = np.array(flines[tables_begin+(table_length*(tab_num-1))+7:tables_begin+(table_length*(tab_num-1))+7+70])
    logTs = [float(s.split()[0]) for s in rows]
    tabvals = [np.asarray(s.split()[1:],dtype=float) for s in rows]
    for i in range(len(tabvals)):
        for j in range(len(tabvals[i])):
            if tabvals[i][j]!=9.999:
                kappas[i][j] = tabvals[i][j]
    if method=="regular":
        interp = RegularGridInterpolator((logTs,logRs), kappas,\
                                     bounds_error=False, fill_value=None)
    else:
        X,Y = np.meshgrid(logTs,logRs)
        tabl = np.vstack((X.flat,Y.flat,kappas.flat)).T
        points = tabl[:,:2]
        values = tabl[:,2]
        print(points)
        print(values)
        points = points[np.isfinite(values)]
        values = values[np.isfinite(values)]
        interp = LinearNDInterpolator(points,values)
    if return_values:
         return interp, logRs,logTs,kappas
    else:
        return interp

def kappa_lookup(T,rho,interp):
    T6 = T/1e6
    logT = np.log10(T)
    R = rho/(T6**3)
    logR = np.log10(R)
    
    return np.power(10,interp((logT,logR)))
