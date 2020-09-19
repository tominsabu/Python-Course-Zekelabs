### HAWT Blade Design
# Development of Blade Element Momentum (BEM) code to calculate the design parameters
# of Horizontal Axis Wind Turbine (HAWT) blades based on input steady wind velocity  
# and required output.  The design parameters also depend on the mechanical and elect-
# -rical efficienty of the turbine components.  

# Importing the required packages
import math    #math contains all math related functions
import pandas as pd    #used to form n-D tables
import numpy as np    #used to form arrays
import matplotlib.pyplot as plt    #pyplot has all the plotting tools
from xfoil import XFoil    #stripped down version of XFoil, available from pypi.org
from xfoil import model    #importing the library 'model' from the package xfoil
from scipy.interpolate import interp1d    #used for 1D numerical interpolation
from scipy.interpolate import InterpolatedUnivariateSpline    #used for interpolation 
import requests
from bs4 import BeautifulSoup
from vpython import*    #vpython is the python graphics package used in this project

# Part1: Scraper
def airfoil_list():
    url="https://m-selig.ae.illinois.edu/ads/coord_database.html"
    response=requests.get(url)    #performs the http request to the website
    # print(response.content)

    soup=BeautifulSoup(response.content,'html.parser')
    # print(soup)

    container=soup.find('div',attrs={'id':'coord_databaseContent'})
    # print(container)
    foils=container.find_all('a', href=True)

    fp=open('airfoil_names.txt','w')
    foil_list=list()
    for foil in foils:
        # foil_list.append(foil.get('href'))
        name=foil.text
        if name.endswith('.dat') or name.endswith('.DAT'):
            fp.write(name) 
            fp.write('\n')
            foil_list.append(name)
        else:
            continue
    fp.close()
    return foil_list
print(airfoil_list())

def naca_list():
    url="https://m-selig.ae.illinois.edu/ads/coord_database.html"
    response=requests.get(url)    #performs the http request to the website
    # print(response.content)

    soup=BeautifulSoup(response.content,'html.parser')
    # print(soup)

    container=soup.find('div',attrs={'id':'coord_databaseContent'})
    # print(container)
    foils=container.find_all('a', href=True)

    fp1=open('nacafoils.txt','w')
    naca_list=list()
    for foil in foils:
        # foil_list.append(foil.get('href'))
        name=foil.text
        if name.startswith('naca') and (name.endswith('.dat') or 
        	name.endswith('.DAT')):
            fp1.write(name)
            fp1.write('\n')
            naca_list.append(name)
        else:
            continue
    fp1.close()
    return naca_list
print(naca_list())

# Downloading the file
def airfoil(airfoil_name):   #  extension also need to be given (.dat or .DAT)
    url="https://m-selig.ae.illinois.edu/ads/coord_database.html"
    response=requests.get(url)    #performs the http request to the website
    # print(response.content)

    soup=BeautifulSoup(response.content,'html.parser')
    # print(soup)

    container=soup.find('div',attrs={'id':'coord_databaseContent'})
    # print(container)
    foils=container.find_all('a', href=True)

    for item in foils:
        name=item.text
        if name==airfoil_name:
            file=requests.get('https://m-selig.ae.illinois.edu/ads/coord/'+name)
            # print(file)
            fp2=open('downloaded_file.txt','wb')
            fp2.write(file.content)
    fp2.close()
    return airfoil_name
print(airfoil('naca4412.dat'))    # when the function is called the airfoil is  
                                  # downloaded to the working directory


# Cleaning the airfoil.txt file for xfoil
stp1=pd.read_fwf('downloaded_file.txt') 
stp1.columns=['y','z']

# last_row={'y':1.000001,'z':0}    
# first_row=pd.DataFrame([pd.Series(last_row)])    
# stp1 = pd.concat([first_row,stp1], ignore_index=True)    # added row [1,0] as 
                                                           # the first row
# stp1= stp1.append(last_row,ignore_index=True)    #added [1,0] as the last row

stp1.at[len(stp1)-1,'z']=stp1.z[0]
print(stp1)

np.savetxt('cleaned_foil.txt', stp1.values, fmt='%1.4e')

# Plotting the airfoil
stp2=pd.read_fwf('cleaned_foil.txt',header=None)
stp2.columns=['y','z']
plt.plot(stp2['y'],stp2['z'])
plt.xlabel('x')
plt.ylabel('y')
plt.show()

# Initial estimate of radius
eta=0.9    # eta is the assumed value for electrical and mechanical efficiences
rho=1.225    # density of air
Cp_initial=0.4    #'Cp_initial' is the expected power coefficient.  As an initial  
                  # guess it is taken as 0.4.  This value will be modified after  
                  # the analysis on the designed blade is carried out.  
u=7    # u is the steady velocity of wind

# 'radius()' will take as input the required power, wind velocity, and  
# power coefficient and returns radius of the blade.  Value of radius, R will  
# be modified after carrying out the analysis.
def radius(P,u,cp):    # P:expected power o/p, U:the wind velocity, cp:expected 
                       # ...power coeff.  
    result=math.sqrt(P/(0.5*cp*eta*rho*3.14*u**3))
    return result

R=radius(10000,u,Cp_initial)    # calls the fuction radius() and assign its value to 
                                # ... variable,R
print('Initial estimate of radius : %.2f meters' %R)



xf1=XFoil()    # invoking the function XFoil() initializes the xfoil object xf1.   
               # A python version of xfoil is available for download with 'pip' and  
               # can be installed by invoking 'pip install xfoil' from the command 
               # prompt.
        

# scraped and cleaned airfoil
airfoil='cleaned_foil.txt'    # airfoil file should be saved in the working directory, 
                              # else the path to the file should be provided


# The function 'foil()' reads the airfoil coordinates from the .txt file and converts  
# the x and y coordinates to numpy arrays.
def foil(airfoil):    # Input is the airfoil.txt file
    foil=pd.read_fwf(airfoil,header=None)  
    foil.columns=['X','Y']    # the columns are given names X and Y respctively
    X=foil['X'].to_numpy()    # to_numpy() converts pandas series X to numpy array
    Y=foil['Y'].to_numpy()
    return X,Y    #returns a list of lists

def unpack_x(airfoil):    # unpacks the list obtained by invoking foil().
    x,_=foil(airfoil)    # '_' is used to store the y cordinates temporarily
    return x  # returns the X coordinates of the airfoil

def unpack_y(airfoil):    # unpack_y() unpacks the list to obtain the y coordinates
    _,y=foil(airfoil)    #  here '_' store the x coordinates temporarily
    return y    # return the Y coordinates of the airfoil

# print(unpack_x('naca0010.txt'))      #x cordinates of the airfoil
# print(unpack_y('naca0010.txt'))      #y cordinates of the airfoil

xf1.airfoil=model.Airfoil(x=unpack_x(airfoil),y=unpack_y(airfoil))   
                # xf1.airfoil tells the xfoil object that an airfoil shape is being 
                # imported.  The model library within the xfoil package contains the 
                # function Airfoil(), which requires the x and y coordinates of the 
                # airfoil.  xf1.airfoil then builds the airfoil shape
  
xf1.repanel(n_nodes=240,cte_ratio=1)    # by default the number of nodes in xfoil 
                                        # would be very less.  nodes can be changed  
                                        # through xf1.repanel(), setting cte_ratio=1 
                                        # will ensure there are more number of nodes
                                        # at leading and trailing edge of airfoil
xf1.Re=1e6    # Reynolds number is set as 1e6
xf1.max_iter=100    #maximum number of iterations set as 100


aoa= np.arange(0,20,0.5)    # create array of angle of attacks in range 0 to 20 
                            # in  steps of 0.5
cl=list()    # creates an empty list of lift coeffs
cd=list()    # creates an empty list of drag coeffs

# for each value in aoa the function xf1.a() calculates along with other qunatities
# the lift and drag coefficients
for value in aoa:
    ans1,ans2,_,_=xf1.a(value)    # lift and drag coefficients.
    cl.append(ans1)    # appends lift coeff to the list cl in each iteration
    cd.append(ans2)    # appends drag coeff to the list cd in each iteration
    
print(cl)
print(cd)


# Tabulating the results obtained from xfoil
# Tabulating aoa,cl,cd obtained from xfoil
data_0={'AoA':aoa,'Cl':cl,'Cd':cd}
df_0=pd.DataFrame(data_0)    # converts data_0 to a pandas dataframe, df_0
df_0['cl/cd']=df_0.Cl/df_0.Cd    # new column 'cl/cd' is created 
df_0.dropna(inplace=True)    # drops any rows which contain 'nan' values.
                             # Some xfoil iterations at higher values of aoa 
                             # might not be converged.  Therefore, they should
                             # be droped to avoid error in further calculations

print(df_0)

# Plotting the airfoil performance curves obtained from xfoil
fig,ax = plt.subplots(1,3)    #subplots(rows, columns)
ax[0].plot(aoa, cl)
ax[0].set_title('Cl Vs AoA')
ax[0].set(xlabel='AoA', ylabel='Cl')

ax[1].plot(aoa, cd)
ax[1].set_title('Cd Vs AoA')
ax[1].set(xlabel='AoA', ylabel='Cd')

ax[2].plot(df_0.AoA,df_0['cl/cd'])
ax[2].set_title('cl/cd Vs AoA')
ax[2].set(xlabel='AoA',ylabel='cl/cd')
plt.show()

# Finding the maximum value of cl/cd
index=df_0[df_0['cl/cd']==df_0['cl/cd'].max()].index[0]   
alpha_opt=df_0.AoA[index]    #alpha correponding to the max value of cl/cd
cl_max=df_0.Cl[index]    #value of lift coeff correspondin to optimum alpha
cd_min=df_0.Cd[index]    #min value of drag coeff correspondin to optimum alpha
print('alpha_opt:',alpha_opt,
	'\t Max lift coeff: %.4f' %cl_max,'\t Min drag coeff: %.4f' %cd_min)


# Dividing the blade into N sections
# When function blade_sections() is called, it divides the blade into 
# the specified number of elements
def blade_sections(N):    # N is the desired number of sections of the blade
    r=[0]    # first radial position is taken as 0 meters from the hub
    element=0
    for i in range(N):
        r.append(element+(R/N))    # radius is divided into N sects, (R/N) gets
        element=element+(R/N)      # appended to the previous value in r
    return r    # returns the list of radaii

N=10    # number of sections taken as 10
r=blade_sections(N)    #calling the function blade_sections()
print('radial positions : ', r)





# Midpoints of blade sections  
# CAUTION: DO NOT CHANGE POSITION OF THIS CELL

# Finding the midpoint of blade sections
r_m=list()    # an empty list is initialised to store the mid section radii
initial=(r[0]+r[1]) / 2    # mid section radii                  
section_width=r[1]-r[0]    # width of blade elements obtained as difference between 
                           # first (index=0) and second (index=1) radial positions
                                       
for count in range(N):               
    r_m.append(initial)
    initial=initial+section_width    # mid sect. radii are found by adding section
                                     # widith to the preceding value
    
print('Mid point of sections : ',r_m)


tip_speed=7    # Tip speed ratio.  A suitable tip speed ratio is selected;  
               # for electrical power generation 4<tip_speed<10.  
B=3    # Number of blades, B is selected as 3

lambda_per_length=tip_speed/R    # tip speed ratio per unit length

# Finding local tip speeds and angles of relative wind speed
lamda=list()    #local tip speeds
phi=list()    #angles of relative wind speed

for count in range(N):
    lamda.append( lambda_per_length*(r_m[count]+(section_width/2)) )    
    # taking r_m here would give the tip speeds at r                                                                 
    phi.append( math.degrees((2/3)*np.arctan(1/lamda[count])) )    
                                                    # phi's are calculated at r
print('Local tip speeds : ',lamda)
print('Relative wind angles : ',phi)



# Finding the local chord distribution
r_new=r[0:]          
r_new.remove(r_new[0])

c=list()
for count in range(N):
    c.append( ((8*3.14*r_new[count])/(B*cl_max))*(1-np.cos(
    	np.radians(phi[count]))) )    
                                    #find the local chord at r itself not r_m
print('List of chords : ',c)

# Taking the pitch angle at the tip as -2
pitch_tip=-2

# Finding the twist angles
gamma=list()
for value in phi:
    gamma.append(value-alpha_opt-pitch_tip)
print('Twist angles : ', gamma)


# Analysing the blade
data0={'r':r_new,'c':c,'gamma':gamma} 
df0=pd.DataFrame(data0)    # The given blade
# with open('table0.tex','w') as tf:    # saves the dataframe as latex table
# 	tf.write(df0.to_latex())
print(df0)

fig1,ax = plt.subplots(1,2)    #subplots(rows, columns)
ax[0].plot(df0.r,df0.c)
ax[0].set_title('chord Vs r')
ax[0].set(xlabel='r', ylabel='chord')

ax[1].plot(df0.r, df0.gamma)
ax[1].set_title('twist Vs r')
ax[1].set(xlabel='r', ylabel='gamma')
plt.show()


#finding the chord lengths corresponding to r_m by interpolation
c_new=c[0:]
c_new.insert(0,c[0])    # chord length at radius = 0 is entered here.  Chord 
                        # length at radius=0  is taken the same as that at the 
                        # succeding radius
f = interp1d(r,c_new)    #inputing x and y into interp1d function
c_m=f(r_m)    # new ys found by inputing new x,c_m's are calculated at each mid-
              # section, starting from the first section.  i.e., the first
              #  section is not disregarded

print(c_m)



data={'r_m':r_m,'c_m':c_m}
df=pd.DataFrame(data)    #table from interpolated data
print(df)


#finding the twist angle, gamma, corresponding to r_m by interpolation
gamma_new=gamma[0:]
gamma_new.insert(0,gamma[0])     #twist at first two positions taken as same

f1 = interp1d(r, gamma_new)    # using interp1d the variation of twist angle, 
                               # gamma_new with radius, r is stored in f1
df['gamma_m']=f1(r_m)    # new values of twist at mid sections, r_m are caculated 
                         #by numerical interpolation
print(df)

# rho=1.225    #density of air
# B=3    #number of blades
# lambda_tip=7    #speed at the tip of the blade
# u=10    #rated wind speed

lamda_new=lamda[0:]
lamda_new.insert(0,0)

f2 = interp1d(r, lamda_new)
df['lambda_m']=f2(r_m)    #tip speed ratios at mid sections got by interpolation

# df['lambda_m']=lambda_tip*(df['r_m']/R)    #not used here

df['sigma_m']=(df['c_m']*B)/(2*3.14*df['r_m'])    #local solidity
df['pitch_m']=df['gamma_m']+pitch_tip    #local pitch
# with open('table1.tex','w') as tf:    #saves pandas dataframe as latex table
# 	tf.write(df.to_latex())
print(df)


# Using InterpolatedUnivariateSpline from scipy.interpolate
order = 1
f_0 = InterpolatedUnivariateSpline(df_0['AoA'],df_0['Cl'], k=order)
# cl_new=f_0(AoA_new)       
# print(cl_new)

f_1 = InterpolatedUnivariateSpline(df_0.AoA,df_0.Cd, k=order)
# cd_new=f_1(AoA_new)       
# print(cd_new)
# print(f_0(0))


# finding the values of 'a_axial','a_radial','phi','angle_attack','Cl',
# 'Cd','F_m','C_T'

a_axial1=0.33    # initial guess for axial induction factor
a_radial1=0    # initial guess for radial induction factor

itervalues=list()    # empty list for storing values of parameters from... 
                     # ... each section of the blade
    
for count in range(N):    # iterate over N sections of the blade
    angle_1=math.degrees(np.arctan( (1-a_axial1)/(df.lambda_m[count]*(
    1+a_radial1)) ))  #relative wind angle
    alpha_m_1=angle_1-df.pitch_m[count]    #angle of attack
  
    Cl=f_0(alpha_m_1) 
    Cd=f_1(alpha_m_1)
   
    f=( (B/2)*(1-(df.r_m[count]/R)) ) / ( (df.r_m[count]/R)*np.sin(math.radians(
    	angle_1)) )
   
    F_m=(2/3.14)*( np.arccos(np.exp(-f)) )    #tip loss correction
    C_T=(df.sigma_m[count]*((1-a_axial1)**2)*(Cl*np.cos(math.radians(angle_1))
        +Cd*np.sin(math.radians(angle_1)))) / (np.sin(math.radians(angle_1))**2)
   
    if C_T< 0.96:
         a_axial_new=1 / (1 + ( ( 4*F_m*(np.sin(math.radians(angle_1)))**2) / (
            df.sigma_m[count]*Cl*np.cos(math.radians(angle_1))) ) )
         a=a_axial1
         b=a_axial_new
         diff=b-a
         a_radial_new=1 / ( ( ( 4*F_m*np.cos(math.radians(angle_1))) / 
            (df.sigma_m[count]*Cl ) ) -1)
         for i in range(8):    # iteration to find the convergence values... 
                                   # ...of a_axial and a_radial
             angle_1=math.degrees(np.arctan( (1-b)/(df.lambda_m[count]*(
             	1+a_radial_new)) ))
             
             alpha_m_1=angle_1-df.pitch_m[count] 
             Cl=f_0(alpha_m_1) 
             Cd=f_1(alpha_m_1)
             f=( (B/2)*(1-(df.r_m[count]/R)) ) / ( (df.r_m[count]/R)*np.sin(
             	math.radians(angle_1)) )
             
             F_m=(2/3.14)*( np.arccos(np.exp(-f)) ) 
             
             C_T=(df.sigma_m[count]*((1-b)**2)*(Cl*np.cos(math.radians(angle_1))
                +Cd*np.sin(math.radians(angle_1)))) / (np.sin(math.radians(
                	angle_1))**2)
             if C_T< 0.96:
                 a=b
                 b=1 / (1 + ( ( 4*F_m*(np.sin(math.radians(angle_1)))**2) /
                  (df.sigma_m[count]*Cl*np.cos(math.radians(angle_1))) ) )
             else:
                 a=b
                 b=(1/F_m)*(0.143+np.sqrt(0.0203-0.6427*(0.889-C_T)))
                
             a_radial_new=1 / ( ( ( 4*F_m*np.cos(math.radians(angle_1))) /
              (df.sigma_m[count]*Cl ) ) -1)
             
             
    else:     #else is executed if C_T > 0.96
        a_axial_new=(1/F_m)*(0.143+np.sqrt(0.0203-0.6427*(0.889-C_T)))
        # print(a_axial_new)
        a=a_axial1
        b=a_axial_new
        diff=b-a
        a_radial_new=1 / ( ( ( 4*F_m*np.cos(math.radians(angle_1))) / 
            (df.sigma_m[count]*Cl ) ) -1)
        
        for i in range(8):
             angle_1=math.degrees(np.arctan( (1-b)/(df.lambda_m[count]*
                (1+a_radial_new)) ))
             alpha_m_1=angle_1-df.pitch_m[count] 
             Cl=f_0(alpha_m_1) 
             Cd=f_1(alpha_m_1)
             f=( (B/2)*(1-(df.r_m[count]/R)) ) / ( (df.r_m[count]/R)*
                np.sin(math.radians(angle_1)) )
             F_m=(2/3.14)*( np.arccos(np.exp(-f)) )
             C_T=(df.sigma_m[count]*(1-b)**2*(Cl*np.cos(math.radians(angle_1))+
                Cd*np.sin(math.radians(angle_1)))) / (np.sin(math.
                radians(angle_1))**2)
             if C_T< 0.96:
                 a=b
                 b=1 / (1 + ( ( 4*F_m*np.sin(math.radians(angle_1))**2) /
                  (df.sigma_m[count]*Cl*np.cos(math.radians(angle_1))) ) )
             else:
                 a=b
                 b=(1/F_m)*(0.143+np.sqrt(0.0203-0.6427*(0.889-C_T)))
             a_radial_new=1 / ( ( ( 4*F_m*np.cos(math.radians(angle_1))) /
              (df.sigma_m[count]*Cl ) ) -1)   
    itervalues.append([b,a_radial_new,angle_1,alpha_m_1,Cl,Cd,F_m,C_T])     
                             # appends outer for loop iteration values
# print(itervalues)
df2=pd.DataFrame(itervalues,columns=['a_axial','a_radial','phi','angle_attack'
    ,'Cl','Cd','F_m','C_T'])
df2.dropna(inplace=True)
print(df2)



df2['lambda_m']=df['lambda_m'] 

#finding the coefficient of power
def coeff_power(tip_speed):
    global Cp
    lambda_tip=tip_speed
    Cp0=8/(lambda_tip*N)
    Cp1=df2.F_m*((np.sin(np.radians(df2.phi)))**2)*(np.cos(np.radians(df2.phi))
        -df2.lambda_m*np.sin(np.radians(df2.phi)))
    Cp2=( np.sin(np.radians(df2.phi)) + df2.lambda_m*np.cos(np.radians(df2.phi)) )*(
     1-(df2.Cd/df2.Cl)*(1/np.tan(np.radians(df2.phi))) )*(df2.lambda_m)**2
    df2['Cp_m']=Cp1*Cp2*Cp0
    Cp=sum(df2.Cp_m)    # coeff of power
    return Cp
print('Coefficient of power = %.4f' %coeff_power(tip_speed))

P=Cp*0.5*rho*3.14*(R**2)*(u**3)    # output power
print('Rated power = %.2f W' %P)

df2_5=df2.drop(['a_axial','a_radial','Cp_m'],axis=1)

# with open('table2.tex','w') as tf:
# 	tf.write(df2_5.to_latex())

df2_6=pd.DataFrame({'axial_induction': df2['a_axial'],'radial_induction'
    :df2['a_radial']})
fig25,ax = plt.subplots(1,2)    #subplots(rows, columns)
ax[0].plot(df.r_m.iloc[1:],df2_6['axial_induction'])
ax[0].set_title('a Vs r')
ax[0].set(xlabel='r', ylabel='a')

ax[1].plot(df.r_m.iloc[1:], df2_6['radial_induction'])
ax[1].set_title(" a'  Vs r ")
ax[1].set(xlabel='r', ylabel="a'")
plt.show()

# with open('table4.tex','w') as tf:
# 	tf.write(df2_6.to_latex())

print(df2_6)
print(df2_5)


# finding lift and drag forces
dr=section_width                 
df3=pd.DataFrame(df2['Cp_m'])
df3['w']=(u*(1-df2.a_axial))/(np.sin(np.radians(df2.phi)))
df3['dFl']=df2.Cl*0.5*rho*df3.w**2*df.c_m*dr      
df3['dFd']=df2.Cd*0.5*rho*df3.w**2*df.c_m*dr
df3['dFn']=df3.dFl*np.cos(np.radians(df2.phi))+df3.dFd*np.sin(
np.radians(df2.phi))     #normal force
df3['dFt']=df3.dFl*np.sin(np.radians(df2.phi))-df3.dFd*np.cos(
np.radians(df2.phi))     #tangential force
df3['dM']=df.r_m*df3.dFt*B 		#torque on turbine at radius r_m
df3['dT']=df3.dFn*B 	#Thrust on blades at a radius r_m
# with open('table3.tex','w') as tf:
# 	tf.write(df3.to_latex())
print(df3)



# Total thrust and torque acting on the the turbine blades
dT_total=sum(df3.dT)
dM_total=sum(df3.dM)
print('Total thrust = %.2f N' %dT_total)           #total thrust
print('Total torque = %.2f Nm' %dM_total)          #total torque



# -----------------------------------------------------------------------------------
# Optimising the tip speed ratio

# Optimizing the tip speed ratio
lambda_list=[entry for entry in range(3,12,1)]
cp_list=list()

for item in lambda_list:
    lambda_per_length=item/R    # tip speed ratio per unit length

    # Finding local tip speeds and angles of relative wind speed
    lamda=list()    #local tip speeds
    phi=list()    #angles of relative wind speed

    for count in range(N):
        lamda.append( lambda_per_length*(r_m[count]+(section_width/2)) )
        phi.append( math.degrees((2/3)*np.arctan(1/lamda[count])) )    

#     print('Local tip speeds : ',lamda)
#     print('Relative wind angles : ',phi)

    # Finding the local chord distribution
    c=list()
    for count in range(N):
        c.append( ((8*3.14*r_m[count])/(B*cl_max))*(1-np.cos(np.radians(
        	phi[count]))) ) 

    # Finding the twist angles
    gamma=list()
    for value in phi:
        gamma.append(value-alpha_opt-pitch_tip)
#     print('Twist angles : ', gamma)

    r_new=r[0:]          
    r_new.remove(r_new[0])

    data0_1={'r':r_new,'c':c,'gamma':gamma} 
    df0_1=pd.DataFrame(data0_1)    

    #finding the chord lengths corresponding to r_m by interpolation
    c_new=c[0:]
    c_new.insert(0,c[0])    # chord length at radius = 0 is entered here.  
                            # Chord length at radius=0 is taken the same as 
                            # that at the succeding radius
    f = interp1d(r,c_new)    #inputing x and y  into interp1d function
    c_m=f(r_m)    #new ys found by inputing new x coordinates

    data={'r_m':r_m,'c_m':c_m}
    df=pd.DataFrame(data)    #table from interpolated data
#     print(df)

    #finding the twist angle, gamma, corresponding to r_m by interpolation
    gamma_new=gamma[0:]
    gamma_new.insert(0,gamma[0])     #twist at first two positions taken as same

    f1 = interp1d(r, gamma_new)    # using interp1d the variation of twist angle, 
                                   # gamma_new with radius, r is stored in f1
    df['gamma_m']=f1(r_m)    # new values of twist at mid sections, r_m are 
                             # caculated by numerical interpolation

    lamda_new=lamda[0:]
    lamda_new.insert(0,0)

    f2 = interp1d(r, lamda_new)
    df['lambda_m']=f2(r_m)    #tip speeds at mid sections got by interpolation

    df['sigma_m']=(df['c_m']*B)/(2*3.14*df['r_m'])    #local solidity
    df['pitch_m']=df['gamma_m']+pitch_tip    #local pitch

    # Using InterpolatedUnivariateSpline from scipy.interpolate
    order = 1
    f_0 = InterpolatedUnivariateSpline(df_0['AoA'],df_0['Cl'], k=order)
    f_1 = InterpolatedUnivariateSpline(df_0.AoA,df_0.Cd, k=order)


    # finding the values of 'a_axial','a_radial','phi','angle_attack','Cl',
    # 'Cd','F_m','C_T'

    a_axial1=0.33    # initial guess for axial induction factor
    a_radial1=0    # initial guess for radial induction factor

    itervalues=list()    # empty list for storing values of parameters from... 
                         # ... each section of the blade

    for count in range(N):    # iterate over N sections of the blade
        angle_1=math.degrees(np.arctan( (1-a_axial1)/(df.lambda_m[count]*(
        1+a_radial1)) ))  #relative wind angle
        alpha_m_1=angle_1-df.pitch_m[count]    #angle of attack

        Cl=f_0(alpha_m_1) 
        Cd=f_1(alpha_m_1)

        f=( (B/2)*(1-(df.r_m[count]/R)) ) / ( (df.r_m[count]/R)*np.sin(math.radians(
        	angle_1)) )

        F_m=(2/3.14)*( np.arccos(np.exp(-f)) )    #tip loss correction
        C_T=(df.sigma_m[count]*((1-a_axial1)**2)*(Cl*np.cos(math.radians(angle_1))
            +Cd*np.sin(math.radians(angle_1)))) / (np.sin(math.radians(angle_1))**2)

        if C_T< 0.96:
             a_axial_new=1 / (1 + ( ( 4*F_m*(np.sin(math.radians(angle_1)))**2) / (
                df.sigma_m[count]*Cl*np.cos(math.radians(angle_1))) ) )
             a=a_axial1
             b=a_axial_new
             diff=b-a
             a_radial_new=1 / ( ( ( 4*F_m*np.cos(math.radians(angle_1))) / 
                (df.sigma_m[count]*Cl ) ) -1)
             for i in range(8):    # iteration to find the convergence values... 
                                       # ...of a_axial and a_radial
                 angle_1=math.degrees(np.arctan( (1-b)/(df.lambda_m[count]*(
                 	1+a_radial_new)) ))

                 alpha_m_1=angle_1-df.pitch_m[count] 
                 Cl=f_0(alpha_m_1) 
                 Cd=f_1(alpha_m_1)
                 f=( (B/2)*(1-(df.r_m[count]/R)) ) / ( (df.r_m[count]/R)*np.sin(
                 	math.radians(angle_1)) )

                 F_m=(2/3.14)*( np.arccos(np.exp(-f)) ) 

                 C_T=(df.sigma_m[count]*((1-b)**2)*(Cl*np.cos(math.radians(angle_1))
                    +Cd*np.sin(math.radians(angle_1)))) / (np.sin(math.radians(
                    	angle_1))**2)
                 if C_T< 0.96:
                     a=b
                     b=1 / (1 + ( ( 4*F_m*(np.sin(math.radians(angle_1)))**2) /
                      (df.sigma_m[count]*Cl*np.cos(math.radians(angle_1))) ) )
                 else:
                     a=b
                     b=(1/F_m)*(0.143+np.sqrt(0.0203-0.6427*(0.889-C_T)))

                 a_radial_new=1 / ( ( ( 4*F_m*np.cos(math.radians(angle_1))) /
                  (df.sigma_m[count]*Cl ) ) -1)


        else:     #else is executed if C_T > 0.96
            a_axial_new=(1/F_m)*(0.143+np.sqrt(0.0203-0.6427*(0.889-C_T)))
            # print(a_axial_new)
            a=a_axial1
            b=a_axial_new
            diff=b-a
            a_radial_new=1 / ( ( ( 4*F_m*np.cos(math.radians(angle_1))) / 
                (df.sigma_m[count]*Cl ) ) -1)

            for i in range(8):
                 angle_1=math.degrees(np.arctan( (1-b)/(df.lambda_m[count]*
                    (1+a_radial_new)) ))
                 alpha_m_1=angle_1-df.pitch_m[count] 
                 Cl=f_0(alpha_m_1) 
                 Cd=f_1(alpha_m_1)
                 f=( (B/2)*(1-(df.r_m[count]/R)) ) / ( (df.r_m[count]/R)*
                    np.sin(math.radians(angle_1)) )
                 F_m=(2/3.14)*( np.arccos(np.exp(-f)) )
                 C_T=(df.sigma_m[count]*(1-b)**2*(Cl*np.cos(math.radians(angle_1))+
                    Cd*np.sin(math.radians(angle_1)))) / (np.sin(math.
                    radians(angle_1))**2)
                 if C_T< 0.96:
                     a=b
                     b=1 / (1 + ( ( 4*F_m*np.sin(math.radians(angle_1))**2) /
                      (df.sigma_m[count]*Cl*np.cos(math.radians(angle_1))) ) )
                 else:
                     a=b
                     b=(1/F_m)*(0.143+np.sqrt(0.0203-0.6427*(0.889-C_T)))
                 a_radial_new=1 / ( ( ( 4*F_m*np.cos(math.radians(angle_1))) /
                  (df.sigma_m[count]*Cl ) ) -1)   
        itervalues.append([b,a_radial_new,angle_1,alpha_m_1,Cl,Cd,F_m,C_T])     
                                 # appends outer for loop iteration values
    # print(itervalues)
    df2=pd.DataFrame(itervalues,columns=['a_axial','a_radial','phi','angle_attack'
        ,'Cl','Cd','F_m','C_T'])
    df2.dropna(inplace=True)
#     print(df2)

    df2['lambda_m']=df['lambda_m'] 
    cp_list.append(coeff_power(item))

cp_list_new=(pd.Series(cp_list)).dropna()
print(cp_list_new)

index2=cp_list_new[cp_list_new==cp_list_new.max()].index[0]
# print(index2)

lambda_opt=lambda_list[index2]
print('Optimum value of tip speed is %.1f' %lambda_opt)
print('Cp at lambda_opt is %.4f' %cp_list_new.max())

fig3,ax = plt.subplots(1,1)    #subplots(rows, columns)
ax.plot(lambda_list,cp_list)
ax.set_title('Cp Vs Lambda')
ax.set(xlabel='Lambda', ylabel='Cp')
plt.show()








# -------------------------------------------------------------------------------
# Linearising the chord


r1=.75*R    # got from literature
r2=.9*R    # got from literature
# print(r1,r2)
f = interp1d(r,c_new)
c1=f(r1)
c2=f(r2)
# print(r1,c1)
# print(r2,c2)

A = np.array([[r1, 1], [r2, 1]])    #[ [x1,y1],[x2,y2] ]
b = np.array([c1, c2])    # other side of equation
a1,a2 = np.linalg.solve(A, b)
print(a1,a2)

# c_lin=a1*r+a2
r_lin=df0.r         #r_lin is same as r
df_lin=pd.DataFrame(r_lin)    
df_lin.columns=['r_lin']
df_lin['c_lin']=a1*r_lin+a2
print(df_lin)


# --------------------------------------------------------------------------------
# Linearising the twist angle


f1 = interp1d(r,gamma_new)
gamma_1=f1(r1)
a3=(gamma_1) / (R-r1)
# gamma_lin=a3*(R-r)    # Linearized twist
df_lin['gamma_lin']=a3*(R-df0.r)
print(df_lin.gamma_lin)

fig1_6,ax = plt.subplots(1,2)    #subplots(rows, columns)
ax[0].plot(df_lin.r_lin,df_lin.c_lin)
ax[0].set_title('chord Vs r')
ax[0].set(xlabel='r', ylabel='chord')

ax[1].plot(df_lin.r_lin,df_lin.gamma_lin)
ax[1].set_title('twist Vs r')
ax[1].set(xlabel='r', ylabel='gamma')
plt.show()






# -------------------------------------------------------------------------------
# Analysing the linearised blade
# --------------------------------------------------------------------------------
tip_speed=lambda_opt    # OPTIMUM TIP SPEED 
B=3    # Number of blades, B is selected as 3
# --------------------------------------------------------------------------------


lambda_per_length=tip_speed/R    # tip speed ratio per unit length

# Finding local tip speeds and angles of relative wind speed
lamda=list()    #local tip speeds
phi=list()    #angles of relative wind speed

for count in range(N):
    lamda.append( lambda_per_length*(r_m[count]+(section_width/2)) )    
                                                                        
                                                                      
    phi.append( math.degrees((2/3)*np.arctan(1/lamda[count])) )    

print('Local tip speeds : ',lamda)
print('Relative wind angles : ',phi)








# Finding the linearized local chord distribution
c=list()
for count in range(N):
    c.append( (a1*r_new[count]) +a2 )    
                                     #find the local chord at r itself not r_m
print('List of chords : ',c)

# Taking the pitch angle at the tip as -2
pitch_tip=-2


gamma=list()
for count in range(N):
    gamma.append(a3*(R-r_new[count]))    #linearized twist
print('Twist angles : ', gamma)


data0_2={'r':r_new,'c':c,'gamma':gamma} 
df0_2=pd.DataFrame(data0_2)    # The given blade
# with open('table0.tex','w') as tf:    # saves the dataframe as latex table
# 	tf.write(df0_2.to_latex())
print(df0_2)

fig1,ax = plt.subplots(1,2)    #subplots(rows, columns)
ax[0].plot(df0_2.r,df0_2.c)
ax[0].set_title('chord Vs r')
ax[0].set(xlabel='r', ylabel='chord')

ax[1].plot(df0_2.r, df0_2.gamma)
ax[1].set_title('twist Vs r')
ax[1].set(xlabel='r', ylabel='gamma')



#finding the chord lengths corresponding to r_m by interpolation
c_new=c[0:]
c_new.insert(0,c[0])    
                        
         
f = interp1d(r,c_new)    
c_m=f(r_m)    
              

print(c_m)


data={'r_m':r_m,'c_m':c_m}
df=pd.DataFrame(data)    #table from interpolated data
print(df)

#finding the twist angle, gamma, corresponding to r_m by interpolation
gamma_new=gamma[0:]
gamma_new.insert(0,gamma[0])     

f1 = interp1d(r, gamma_new)    
df['gamma_m']=f1(r_m)    
print(df)


# rho=1.225    #density of air
# B=3    #number of blades
# lambda_tip=7    #speed at the tip of the blade
# u=10    #rated wind speed

lamda_new=lamda[0:]
lamda_new.insert(0,0)

f2 = interp1d(r, lamda_new)
df['lambda_m']=f2(r_m)    #tip speed ratios at mid sections got by interpolation

# df['lambda_m']=lambda_tip*(df['r_m']/R)    #not used here

df['sigma_m']=(df['c_m']*B)/(2*3.14*df['r_m'])    #local solidity
df['pitch_m']=df['gamma_m']+pitch_tip    #local pitch
# with open('table1.tex','w') as tf:    #saves pandas dataframe as a latex table
# 	tf.write(df.to_latex())
print(df)

# Using InterpolatedUnivariateSpline from scipy.interpolate
order = 1
f_0 = InterpolatedUnivariateSpline(df_0['AoA'],df_0['Cl'], k=order)
# cl_new=f_0(AoA_new)       
# print(cl_new)

f_1 = InterpolatedUnivariateSpline(df_0.AoA,df_0.Cd, k=order)
# cd_new=f_1(AoA_new)       
# print(cd_new)
# print(f_0(0))




# finding the values of 'a_axial','a_radial','phi','angle_attack','Cl',
# 'Cd','F_m','C_T'

a_axial1=0.33    # initial guess for axial induction factor
a_radial1=0    # initial guess for radial induction factor

itervalues=list()    # empty list for storing values of parameters from... 
                     # ... each section of the blade
    
for count in range(N):    # iterate over N sections of the blade
    angle_1=math.degrees(np.arctan( (1-a_axial1)/(df.lambda_m[count]*(
    1+a_radial1)) ))  #relative wind angle
    alpha_m_1=angle_1-df.pitch_m[count]    #angle of attack
  
    Cl=f_0(alpha_m_1) 
    Cd=f_1(alpha_m_1)
   
    f=( (B/2)*(1-(df.r_m[count]/R)) ) / ( (df.r_m[count]/R)*np.sin(
    	math.radians(angle_1)) )
   
    F_m=(2/3.14)*( np.arccos(np.exp(-f)) )    #tip loss correction
    C_T=(df.sigma_m[count]*((1-a_axial1)**2)*(Cl*np.cos(math.radians(angle_1))
        +Cd*np.sin(math.radians(angle_1)))) / (np.sin(math.radians(angle_1))**2)
   
    if C_T< 0.96:
         a_axial_new=1 / (1 + ( ( 4*F_m*(np.sin(math.radians(angle_1)))**2) / (
            df.sigma_m[count]*Cl*np.cos(math.radians(angle_1))) ) )
         a=a_axial1
         b=a_axial_new
         diff=b-a
         a_radial_new=1 / ( ( ( 4*F_m*np.cos(math.radians(angle_1))) / 
            (df.sigma_m[count]*Cl ) ) -1)
         for i in range(8):    # iteration to find the convergence values... 
                                   # ...of a_axial and a_radial
             angle_1=math.degrees(np.arctan( (1-b)/(df.lambda_m[count]*(
             	1+a_radial_new)) ))
             
             alpha_m_1=angle_1-df.pitch_m[count] 
             Cl=f_0(alpha_m_1) 
             Cd=f_1(alpha_m_1)
             f=( (B/2)*(1-(df.r_m[count]/R)) ) / ( (df.r_m[count]/R)*np.sin(
             	math.radians(angle_1)) )
             
             F_m=(2/3.14)*( np.arccos(np.exp(-f)) ) 
             
             C_T=(df.sigma_m[count]*((1-b)**2)*(Cl*np.cos(math.radians(angle_1))
                +Cd*np.sin(math.radians(angle_1)))) / (np.sin(math.radians(
                	angle_1))**2)
             if C_T< 0.96:
                 a=b
                 b=1 / (1 + ( ( 4*F_m*(np.sin(math.radians(angle_1)))**2) /
                  (df.sigma_m[count]*Cl*np.cos(math.radians(angle_1))) ) )
             else:
                 a=b
                 b=(1/F_m)*(0.143+np.sqrt(0.0203-0.6427*(0.889-C_T)))
                
             a_radial_new=1 / ( ( ( 4*F_m*np.cos(math.radians(angle_1))) /
              (df.sigma_m[count]*Cl ) ) -1)
             
             
    else:     #else is executed if C_T > 0.96
        a_axial_new=(1/F_m)*(0.143+np.sqrt(0.0203-0.6427*(0.889-C_T)))
        # print(a_axial_new)
        a=a_axial1
        b=a_axial_new
        diff=b-a
        a_radial_new=1 / ( ( ( 4*F_m*np.cos(math.radians(angle_1))) / 
            (df.sigma_m[count]*Cl ) ) -1)
        
        for i in range(8):
             angle_1=math.degrees(np.arctan( (1-b)/(df.lambda_m[count]*
                (1+a_radial_new)) ))
             alpha_m_1=angle_1-df.pitch_m[count] 
             Cl=f_0(alpha_m_1) 
             Cd=f_1(alpha_m_1)
             f=( (B/2)*(1-(df.r_m[count]/R)) ) / ( (df.r_m[count]/R)*
                np.sin(math.radians(angle_1)) )
             F_m=(2/3.14)*( np.arccos(np.exp(-f)) )
             C_T=(df.sigma_m[count]*(1-b)**2*(Cl*np.cos(math.radians(angle_1))+
                Cd*np.sin(math.radians(angle_1)))) / (np.sin(math.
                radians(angle_1))**2)
             if C_T< 0.96:
                 a=b
                 b=1 / (1 + ( ( 4*F_m*np.sin(math.radians(angle_1))**2) /
                  (df.sigma_m[count]*Cl*np.cos(math.radians(angle_1))) ) )
             else:
                 a=b
                 b=(1/F_m)*(0.143+np.sqrt(0.0203-0.6427*(0.889-C_T)))
             a_radial_new=1 / ( ( ( 4*F_m*np.cos(math.radians(angle_1))) /
              (df.sigma_m[count]*Cl ) ) -1)   
    itervalues.append([b,a_radial_new,angle_1,alpha_m_1,Cl,Cd,F_m,C_T])     
                             # appends outer for loop iteration values
# print(itervalues)
df2=pd.DataFrame(itervalues,columns=['a_axial','a_radial','phi','angle_attack'
    ,'Cl','Cd','F_m','C_T'])
df2.dropna(inplace=True)
print(df2)



df2['lambda_m']=df['lambda_m'] 

#finding the coefficient of power
def coeff_power(tip_speed):
    global Cp
    lambda_tip=tip_speed
    Cp0=8/(lambda_tip*N)
    Cp1=df2.F_m*((np.sin(np.radians(df2.phi)))**2)*(np.cos(np.radians(df2.phi))
        -df2.lambda_m*np.sin(np.radians(df2.phi)))
    Cp2=( np.sin(np.radians(df2.phi)) + df2.lambda_m*np.cos(np.radians(df2.phi)) )*(
     1-(df2.Cd/df2.Cl)*(1/np.tan(np.radians(df2.phi))) )*(df2.lambda_m)**2
    df2['Cp_m']=Cp1*Cp2*Cp0
    Cp=sum(df2.Cp_m)    # coeff of power
    return Cp
print('Coefficient of power = %.4f' %coeff_power(tip_speed))

P=Cp*0.5*rho*3.14*(R**2)*(u**3)    # output power
print('Rated power = %.2f W' %P)

df2_5=df2.drop(['a_axial','a_radial','Cp_m'],axis=1)
# with open('table2.tex','w') as tf:
# 	tf.write(df2_5.to_latex())

df2_6=pd.DataFrame({'axial_induction': df2['a_axial'],'radial_induction'
    :df2['a_radial']})
# with open('table4.tex','w') as tf:
# 	tf.write(df2_6.to_latex())

fig26,ax = plt.subplots(1,2)    #subplots(rows, columns)
ax[0].plot(df_lin.r_lin,df2_6['axial_induction'])
ax[0].set_title('a Vs r')
ax[0].set(xlabel='r', ylabel='a')

ax[1].plot(df_lin.r_lin, df2_6['radial_induction'])
ax[1].set_title(" a'  Vs r ")
ax[1].set(xlabel='r', ylabel="a'")
plt.show()

print(df2_6)
print(df2_5)









# --------------------------------------------------------------------------------

# cp-lambda performance curve, showing variation with solidity




lambda_list=[item for item in range(3,20,1)]
blades=[1,2,3,4,5]

cp_list_list=list()

ct_list_list=list()
for B in blades:
    cp_list=list()
    
    ct_list=list()
    for tip_speed in lambda_list:

        lambda_per_length=tip_speed/R    # tip speed ratio per unit length

        # Finding local tip speeds and angles of relative wind speed
        lamda=list()    #local tip speeds
        phi=list()    #angles of relative wind speed

        for count in range(N):
            lamda.append( lambda_per_length*(r_m[count]+(section_width/2)) )    
                                                                                
            phi.append( math.degrees((2/3)*np.arctan(1/lamda[count])) )     



        # Finding the linearized local chord distribution
        c=list()
        for count in range(N):
            c.append( (a1*r_new[count]) +a2 )    


        # Taking the pitch angle at the tip as -2
        pitch_tip=-2

        # Use the linearized twist angle here
        gamma=list()
        for count in range(N):
            gamma.append(a3*(R-r_new[count]))    #linearized twist



        #finding the chord lengths corresponding to r_m by interpolation
        c_new=c[0:]
        c_new.insert(0,c[0])    

        f = interp1d(r,c_new)    
        c_m=f(r_m)   

#         print(c_m)


        data={'r_m':r_m,'c_m':c_m}
        df=pd.DataFrame(data)    #table from interpolated data
#         print(df)

        #finding the twist angle, gamma, corresponding to r_m by interpolation
        gamma_new=gamma[0:]
        gamma_new.insert(0,gamma[0])     

        f1 = interp1d(r, gamma_new)    
        df['gamma_m']=f1(r_m)   


        lamda_new=lamda[0:]
        lamda_new.insert(0,0)

        f2 = interp1d(r, lamda_new)
        df['lambda_m']=f2(r_m)    

        # df['lambda_m']=lambda_tip*(df['r_m']/R)    #not used here

        df['sigma_m']=(df['c_m']*B)/(2*3.14*df['r_m'])    #local solidity
        df['pitch_m']=df['gamma_m']+pitch_tip    #local pitch
       

        # Using InterpolatedUnivariateSpline from scipy.interpolate
        order = 1
        f_0 = InterpolatedUnivariateSpline(df_0['AoA'],df_0['Cl'], k=order)
        # cl_new=f_0(AoA_new)       
        # print(cl_new)

        f_1 = InterpolatedUnivariateSpline(df_0.AoA,df_0.Cd, k=order)
        # cd_new=f_1(AoA_new)       
        # print(cd_new)
        # print(f_0(0))




        a_axial1=0.33    # initial guess for axial induction factor
        a_radial1=0    # initial guess for radial induction factor

        itervalues=list()    

        for count in range(N):    # iterate over N sections of the blade
            angle_1=math.degrees(np.arctan( (1-a_axial1)/(df.lambda_m[count]*(
            1+a_radial1)) ))  #relative wind angle
            alpha_m_1=angle_1-df.pitch_m[count]    #angle of attack

            Cl=f_0(alpha_m_1) 
            Cd=f_1(alpha_m_1)

            f=( (B/2)*(1-(df.r_m[count]/R)) ) / ( (df.r_m[count]/R)*np.sin(
            	math.radians(angle_1)) )

            F_m=(2/3.14)*( np.arccos(np.exp(-f)) )    #tip loss correction
            C_T=(df.sigma_m[count]*((1-a_axial1)**2)*(Cl*np.cos(math.radians(angle_1))
                +Cd*np.sin(math.radians(angle_1)))) / (np.sin(math.radians(angle_1))**2)
            
            if C_T< 0.9600:
                 a_axial_new=1 / (1 + ( ( 4*F_m*(np.sin(math.radians(angle_1)))**2) / (
                    df.sigma_m[count]*Cl*np.cos(math.radians(angle_1))) ) )
                 a=a_axial1
                 b=a_axial_new
                 diff=b-a
                 a_radial_new=1 / ( ( ( 4*F_m*np.cos(math.radians(angle_1))) / 
                    (df.sigma_m[count]*Cl ) ) -1)
                 for i in range(8):    # iteration to find the convergence values... 
                                           # ...of a_axial and a_radial
                     angle_1=math.degrees(np.arctan( (1-b)/(df.lambda_m[count]*(
                     	1+a_radial_new)) ))

                     alpha_m_1=angle_1-df.pitch_m[count] 
                     Cl=f_0(alpha_m_1) 
                     Cd=f_1(alpha_m_1)
                     f=( (B/2)*(1-(df.r_m[count]/R)) ) / ( (df.r_m[count]/R)*np.sin(
                     	math.radians(angle_1)) )

                     F_m=(2/3.14)*( np.arccos(np.exp(-f)) ) 

                     C_T=(df.sigma_m[count]*((1-b)**2)*(Cl*np.cos(math.radians(angle_1))
                        +Cd*np.sin(math.radians(angle_1)))) / (np.sin(math.radians(
                        	angle_1))**2)
                     if C_T< 0.960:
                         a=b
                         b=1 / (1 + ( ( 4*F_m*(np.sin(math.radians(angle_1)))**2) /
                          (df.sigma_m[count]*Cl*np.cos(math.radians(angle_1))) ) )
                     elif C_T>0.960:
                         a=b
                         b=(1/F_m)*(0.143+np.sqrt(0.0203-0.6427*(0.889-C_T)))

                     a_radial_new=1 / ( ( ( 4*F_m*np.cos(math.radians(angle_1))) /
                      (df.sigma_m[count]*Cl ) ) -1)


            elif C_T>0.960:     #else is executed if C_T > 0.96
                a_axial_new=(1/F_m)*(0.143+np.sqrt(0.0203-0.6427*(0.889-C_T)))
                # print(a_axial_new)
                a=a_axial1
                b=a_axial_new
                diff=b-a
                a_radial_new=1 / ( ( ( 4*F_m*np.cos(math.radians(angle_1))) / 
                    (df.sigma_m[count]*Cl ) ) -1)

                for i in range(8):
                     angle_1=math.degrees(np.arctan( (1-b)/(df.lambda_m[count]*
                        (1+a_radial_new)) ))
                     alpha_m_1=angle_1-df.pitch_m[count] 
                     Cl=f_0(alpha_m_1) 
                     Cd=f_1(alpha_m_1)
                     f=( (B/2)*(1-(df.r_m[count]/R)) ) / ( (df.r_m[count]/R)*
                        np.sin(math.radians(angle_1)) )
                     F_m=(2/3.14)*( np.arccos(np.exp(-f)) )
                     C_T=(df.sigma_m[count]*(1-b)**2*(Cl*np.cos(math.radians(angle_1))+
                        Cd*np.sin(math.radians(angle_1)))) / (np.sin(math.
                        radians(angle_1))**2)
                     if C_T< 0.96:
                         a=b
                         b=1 / (1 + ( ( 4*F_m*np.sin(math.radians(angle_1))**2) /
                          (df.sigma_m[count]*Cl*np.cos(math.radians(angle_1))) ) )
                     elif C_T>0.960:
                         a=b
                         b=(1/F_m)*(0.143+np.sqrt(0.0203-0.6427*(0.889-C_T)))
                     a_radial_new=1 / ( ( ( 4*F_m*np.cos(math.radians(angle_1))) /
                      (df.sigma_m[count]*Cl ) ) -1)   
            itervalues.append([b,a_radial_new,angle_1,alpha_m_1,Cl,Cd,F_m,C_T])     
                          # appends outer for loop iteration values
        # print(itervalues)
        df2=pd.DataFrame(itervalues,columns=['a_axial','a_radial','phi','angle_attack'
            ,'Cl','Cd','F_m','C_T'])
        df2.dropna(inplace=True)
#         print(df2)
        


        df2['lambda_m']=df['lambda_m'] 
       
        cp_list.append(coeff_power(tip_speed))

        
    cp_list_list.append(cp_list)
    ct_list_list.append(ct_list)
    
            
labels=['B=1','B=2','B=3','B=4','B=5']

for item,label in zip(cp_list_list,labels):
    plt.plot(lambda_list, item, label=label)
plt.xlabel('Tip speed ratio')
plt.ylabel('Coefficient of Power')
plt.legend()
plt.show()






# --------------------------------------------------------------------------------
# Graphics


airfoil=pd.read_fwf('cleaned_foil.txt',header=None)     #best way to read text files
airfoil.columns=['X','Y']
shape=airfoil[['X','Y']].values.tolist()       

                
finer=1               #divides each section of the blade into finer subsections. 
                              # Greater the value of finer, the smoother the blade.

s1=section_width/finer
paths=list()

for x in range(finer):
    pth=[vector((s1+(x-1)*s1),0,0),vector((s1+x*s1),0,0)]
    paths.append(pth)

# p1=[vector(0,0,0),vector(s1,0,0)]
# p2=[vector(s1,0,0),vector(s1+s1,0,0)]

sect=list()
for count0 in range(len(df0)):
    ext=list()

    gamma=df0.gamma[count0]
    gamma_x=list()

    if count0==0:
        gamma_x_1=0
        for i in range(finer):
            gamma_x_2=gamma_x_1+gamma/finer
            gamma_x.append(gamma_x_2)
            gamma_x_1=gamma_x_2
    else:
        for i in range(finer):

            gamma_x_2=gamma_x_1-gamma/finer
            gamma_x.append(gamma_x_2)
            gamma_x_1=gamma_x_2

    count=0
    for value in paths:
        ext_x=extrusion(shape=shape,path=value)
        ext_x.rotate(angle=np.radians(-gamma_x[count]))
        if count==(finer-1):
            if count0<9:
                ext_x.width=(df0.c[count0]+df0.c[count0+1])/2
            elif count0==9:
                ext_x.width=df0.c[count0]
        else:
            ext_x.width=df0.c[count0]
        # ext_x.length=0.0649
        ext.append(ext_x)
        count+=1

    sect1 = compound(ext,pos=vector(df0.r[count0]-3.5,1,3))
    blade=sect.append(sect1)


# Again doing with finer = 15

finer=15               

s1=section_width/finer
paths=list()

for x in range(finer):
    pth=[vector((s1+(x-1)*s1),0,0),vector((s1+x*s1),0,0)]
    paths.append(pth)

# p1=[vector(0,0,0),vector(s1,0,0)]
# p2=[vector(s1,0,0),vector(s1+s1,0,0)]

sect=list()
for count0 in range(len(df0)):
    ext=list()

    gamma=df0.gamma[count0]
    gamma_x=list()

    if count0==0:
        gamma_x_1=0
        for i in range(finer):
            gamma_x_2=gamma_x_1+gamma/finer
            gamma_x.append(gamma_x_2)
            gamma_x_1=gamma_x_2
    else:
        for i in range(finer):

            gamma_x_2=gamma_x_1-gamma/finer
            gamma_x.append(gamma_x_2)
            gamma_x_1=gamma_x_2

    count=0
    for value in paths:
        ext_x=extrusion(shape=shape,path=value)
        ext_x.rotate(angle=np.radians(-gamma_x[count]))
        if count==(finer-1):
            if count0<9:
                ext_x.width=(df0.c[count0]+df0.c[count0+1])/2
            elif count0==9:
                ext_x.width=df0.c[count0]
        else:
            ext_x.width=df0.c[count0]
        # ext_x.length=0.0649
        ext.append(ext_x)
        count+=1

    sect1 = compound(ext,pos=vector(df0.r[count0]-3.5,0,3))
    blade=sect.append(sect1)



# Linearized blade

   
finer=1               

s1=section_width/finer
paths=list()

for x in range(finer):
    pth=[vector((s1+(x-1)*s1),0,0),vector((s1+x*s1),0,0)]
    paths.append(pth)

# p1=[vector(0,0,0),vector(s1,0,0)]
# p2=[vector(s1,0,0),vector(s1+s1,0,0)]

sect=list()
for count0 in range(len(df0_2)):
    ext=list()

    gamma=df0_2.gamma[count0]
    gamma_x=list()

    if count0==0:
        gamma_x_1=0
        for i in range(finer):
            gamma_x_2=gamma_x_1+gamma/finer
            gamma_x.append(gamma_x_2)
            gamma_x_1=gamma_x_2
    else:
        for i in range(finer):

            gamma_x_2=gamma_x_1-gamma/finer
            gamma_x.append(gamma_x_2)
            gamma_x_1=gamma_x_2

    count=0
    for value in paths:
        ext_x=extrusion(shape=shape,path=value)
        ext_x.rotate(angle=np.radians(-gamma_x[count]))
        if count==(finer-1):
            if count0<9:
                ext_x.width=(df0_2.c[count0]+df0_2.c[count0+1])/2
            elif count0==9:
                ext_x.width=df0_2.c[count0]
        else:
            ext_x.width=df0_2.c[count0]
        # ext_x.length=0.0649
        ext.append(ext_x)
        count+=1

    sect1 = compound(ext,pos=vector(df0_2.r[count0]-3.5,-1,3))
    blade=sect.append(sect1)




finer=15               

s1=section_width/finer
paths=list()

for x in range(finer):
    pth=[vector((s1+(x-1)*s1),0,0),vector((s1+x*s1),0,0)]
    paths.append(pth)

# p1=[vector(0,0,0),vector(s1,0,0)]
# p2=[vector(s1,0,0),vector(s1+s1,0,0)]

sect=list()
for count0 in range(len(df0_2)):
    ext=list()

    gamma=df0_2.gamma[count0]
    gamma_x=list()

    if count0==0:
        gamma_x_1=0
        for i in range(finer):
            gamma_x_2=gamma_x_1+gamma/finer
            gamma_x.append(gamma_x_2)
            gamma_x_1=gamma_x_2
    else:
        for i in range(finer):

            gamma_x_2=gamma_x_1-gamma/finer
            gamma_x.append(gamma_x_2)
            gamma_x_1=gamma_x_2

    count=0
    for value in paths:
        ext_x=extrusion(shape=shape,path=value)
        ext_x.rotate(angle=np.radians(-gamma_x[count]))
        if count==(finer-1):
            if count0<9:
                ext_x.width=(df0_2.c[count0]+df0_2.c[count0+1])/2
            elif count0==9:
                ext_x.width=df0_2.c[count0]
        else:
            ext_x.width=df0_2.c[count0]
        # ext_x.length=0.0649
        ext.append(ext_x)
        count+=1

    sect1 = compound(ext,pos=vector(df0_2.r[count0]-3.5,-2,3))
    blade=sect.append(sect1)
