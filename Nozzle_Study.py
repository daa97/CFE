from Nozzle_Main import *

#default_q_regen = get_q_nozzle() # default_q_regen is 794 W
#ic(default_q_regen)
#################################### TRADE STUDIES: WHAT HAPPENS IF I CHANGE ______? ################################

#Optimization

chan_geo1 = {"k_w": 398, "t_w": .0254, "chan_width": .0254 * 1, "chan_depth": .0254 * .01, "chan_gap":.0254 * .01} #fat and short #previously 4030 before wall thickness
#chan_geo12 = {"k_w": 398, "t_w": .00254, "chan_width": .0254 * 1, "chan_depth": .0254 * .01, "chan_gap":.0254 * .01} #fat and short 4090
#chan_geo13 = {"k_w": 398, "t_w": .00254, "chan_width": .0254 * 1, "chan_depth": .0254 * 1/8, "chan_gap":.0254 * .01} #  fat fat skinny wall #4088
#chan_geo14 = {"k_w": 398, "t_w": .00254, "chan_width": .0254 * 1/2, "chan_depth": .0254 * 1/2, "chan_gap":.0254 * .01} #  fat fat skinny wall #4598
#chan_geo2 = {"k_w": 398, "t_w": .0254, "chan_width": .0254 * .01, "chan_depth": .0254 * 1, "chan_gap":.0254 * .01} #skinny and tall # 2540 before wall thickness
#chan_geo22 = {"k_w": 398, "t_w": .00254, "chan_width": .0254 * .01, "chan_depth": .0254 * 1, "chan_gap":.0254 * .01} #skinny and tall # 2577

q_regen = get_q_nozzle(chan_geo=chan_geo1, e=200)
ic(q_regen)

'''chan_geo3 = {"k_w": 398, "t_w": .00254, "chan_width": .0254 * 1/8, "chan_depth": .0254 * 1/8, "chan_gap":.0254 * .01} # medium medium skinny wall # 4730
chan_geo32 = {"k_w": 398, "t_w": .00254, "chan_width": .0254 * 1, "chan_depth": .0254 * 1, "chan_gap":.0254 * .01} #  fat fat skinny wall 4078
chan_geo33 = {"k_w": 398, "t_w": .00254, "chan_width": .0254 * 3/4, "chan_depth": .0254 * 3/4, "chan_gap":.0254 * .01} # medium medium skinny wall #4593
chan_geo34 = {"k_w": 398, "t_w": .00254, "chan_width": .0254 * 1/16, "chan_depth": .0254 * 1/16, "chan_gap":.0254 * .01} # medium medium skinny wall #



min_value = 1/16
max_value = 3/4
data_pts = 10

heat = []

for a in np.linspace(min_value, max_value, data_pts):
    # Change parameter
    # Compute New Nozzle Geometry
    new_chan_geo = {"k_w": 398, "t_w": .00254, "chan_width": a*.0254, "chan_depth": a*.0254, "chan_gap": .0254*.01}
    q_regen = get_q_nozzle(chan_geo=new_chan_geo, e=200)
    ic(q_regen)
    # Assignment
    heat = np.append(heat, q_regen)

# PLOTS
plt.plot((np.linspace(min_value, max_value, data_pts)), heat)
plt.xlabel("Nozzle Wall Thermal Conducitivty (W/m^2K")
plt.ylabel("Total Nozzle Heat Transfer (W)")
plt.show()
ic(heat)'''

'''array([4410.60212633, 4687.16782936, 4842.782945  , 4771.01172645,
                 4890.79120033, 4996.20718396, 4789.2711922 , 4880.25595924,
                 4815.0808836 , 4593.88422793])'''


#================================== STANDARD 1-D ISENTROPIC VERSUS ""+ NASA CEA =========================
'''Nozzle.noz_heat.oneDIsentropic(Nozzle,showplots=True)
plt.show()'''

#============================================ NUMBER OF DATA POINTS  ================================================
'''min_value = 10
max_value = 40
data_pts = 4

heat_N = []

for a in np.linspace(min_value, max_value, data_pts):
    # Change parameter
    N_data_pts = np.floor(a)
    # Compute New Nozzle Geometry
    q_regen = get_q_nozzle(data_pts=N_data_pts)
    #Assignment
    heat_N = np.append(heat_N, q_regen)
    
#PLOTS
plt.plot((np.linspace(min_value, max_value, data_pts)), heat_N)
plt.xlabel("Number of Data points")
plt.ylabel("Total Nozzle Heat Transfer (W)")
plt.show()'''
'''#============================================ AREA RATIO  ================================================
min_value = 50
max_value = 300
data_pts = 10

heat_e = []

for a in np.linspace(min_value, max_value, data_pts):
    # Change parameter
    # Compute New Nozzle Geometry
    q_regen = get_q_nozzle(e=a)
    #Assignment
    heat_e = np.append(heat_e, q_regen)
    
#PLOTS
plt.plot((np.linspace(min_value, max_value, data_pts)), heat_e)
plt.xlabel("Nozzle Area Ratio")
plt.ylabel("Total Nozzle Heat Transfer (W)")
plt.show()'''


'''#============================================ THERMAL CONDUCTIVITY ================================================
# Zirconium Carbide has superior resistance to fissile products at high temperatures. It's k= 20 W/mk
min_value = 10
max_value = 400
data_pts = 10

heat_k_w=[]

for a in np.linspace(min_value, max_value, data_pts):
    # Change parameter
    # Compute New Nozzle Geometry
    new_chan_geo = {"k_w": a, "t_w": .00254, "chan_width": default_w, "chan_depth": default_h, "chan_gap":default_g}
    q_regen = get_q_nozzle(chan_geo=new_chan_geo)
    #Assignment
    heat_k_w = np.append(heat_k_w, q_regen)
    
#PLOTS
plt.plot((np.linspace(min_value, max_value, data_pts)), heat_k_w)
plt.xlabel("Nozzle Wall Thermal Conducitivty (W/m^2K")
plt.ylabel("Total Nozzle Heat Transfer (W)")
plt.show(block=False)'''



#============================================ WALL THICKNESS ================================================


'''min_value =.0254*.01
max_value = .0254 # between ten thou and 1 inch
data_pts = 10

heat_t_w = []

for a in np.linspace(min_value, max_value, data_pts):
    # Change parameter
    # Compute New Nozzle Geometry
    new_chan_geo = {"k_w": 398, "t_w": .00254, "chan_width": default_w, "chan_depth": default_h, "chan_gap": default_g}
    q_regen = get_q_nozzle(chan_geo=new_chan_geo)
    # Assignment
    heat_t_w = np.append(heat_t_w, q_regen)

# PLOTS
plt.plot((np.linspace(min_value, max_value, data_pts)), heat_t_w)
plt.xlabel("Wall Thickness (m)")
plt.ylabel("Total Nozzle Heat Transfer (W)")
plt.show(block=False)
plt.savefig('t_w.png')
'''

'''#============================================ COOLANT VELOCITY  ================================================
# (by changing flow rate channel area) from ten thou by ten thou to 1/8" by 1/8" (changing channel dimensions)


min_value =.0254*.01
max_value = .0254*.25 # between 10 thou and a quarter-inch
data_pts = 10

heat_chan_width = []

for a in np.linspace(min_value, max_value, data_pts):
    # Change parameter
    # Compute New Nozzle Geometry
    new_chan_geo = {"k_w": 398, "t_w": .00254, "chan_width": a, "chan_depth": default_h, "chan_gap": default_g}
    q_regen = get_q_nozzle(chan_geo=new_chan_geo)
    # Assignment
    ic(q_regen)
    heat_chan_width = np.append(heat_chan_width, q_regen)

# PLOTS
plt.plot((np.linspace(min_value, max_value, data_pts)), heat_chan_width, '-ko')
plt.xlabel("Channel Width (m)")
plt.ylabel("Total Nozzle Heat Transfer (W)")
plt.show()
plt.savefig('c_w.png')'''

#====================================================== Channel Height=============================================

'''min_value =.0254*.01
max_value = .0254 # between 10 thou and an inch
data_pts = 10

heat_chan_h = []

for a in np.linspace(min_value, max_value, data_pts):
    # Change parameter
    # Compute New Nozzle Geometry
    new_chan_geo = {"k_w": 398, "t_w": .00254, "chan_width": default_w, "chan_depth": a, "chan_gap": default_g}
    q_regen = get_q_nozzle(chan_geo=new_chan_geo)
    # Assignment
    heat_chan_h = np.append(heat_chan_h, q_regen)

# PLOTS
plt.figure()
plt.plot((np.linspace(min_value, max_value, data_pts)), heat_chan_h, '-ko')
plt.xlabel("Channel Height (m)")
plt.ylabel("Total Nozzle Heat Transfer (W)")
plt.show()

#====================================================== Channel Gap =================================================

min_value =.0254*.01
max_value = .0254 # between 10 thou and an inch
data_pts = 10

heat_chan_g = []

for a in np.linspace(min_value, max_value, data_pts):
    # Change parameter
    # Compute New Nozzle Geometry
    new_chan_geo = {"k_w": 398, "t_w": .00254, "chan_width": default_w, "chan_depth": default_h, "chan_gap": a}
    q_regen = get_q_nozzle(chan_geo=new_chan_geo)
    # Assignment
    heat_chan_g = np.append(heat_chan_g, q_regen)

# PLOTS
plt.figure()
plt.plot((np.linspace(min_value, max_value, data_pts)), heat_chan_g,'-bo')
plt.xlabel("Channel Gap (m)")
plt.ylabel("Total Nozzle Heat Transfer (W)")
plt.savefig('c_g.png')
plt.show()'''
