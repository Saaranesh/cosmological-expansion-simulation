#Importing libraries
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import numpy as np
from scipy.integrate import quad

def friedmann(z_prime, Or, Om0, Ok0, Ode0): #Establishes and returns integrand of hubble parameter
    zp1 = z_prime + 1
    return 1/(zp1 * (np.sqrt((Or * zp1**4) + (Om0 * zp1**3) + (Ok0 * zp1**2) + Ode0)))

  def age(z, H0, Om0, Ode0): #Calculates time elapsed since big bang at a point in time with redshift z
    Or = 8.24e-5 #Estimate from Planck satellite
    Om0 = Om0 - Or/2
    Ode0 = Ode0 - Or/2
    Ok0 = 1-(Om0+Ode0+Or) #This should result in 0, as we assume a spatially flat universe (k = 0)
    C = 978 #Conversion constant between 1/H0 and time in Gyr
    integrate_result, error_estimate = quad(friedmann, z, np.inf, args=(Or, Om0, Ok0, Ode0)) #Integral evaluator executes computation of integral in hubble parameter with limits infinity (big bang) and z
    t = np.abs((C/H0) * integrate_result) #Computes t in gigayears
    return t

def z_at_time(t, minz, maxz, H0, Om0, Ode0): #Calculates value of cosmological redshift at a point in time, t
    if (age(minz, H0, Om0, Ode0) - t) * (age(maxz, H0, Om0, Ode0) - t) > 0: #This will be negative if the target z value lies in the range of minz to maxz
        return None #hence validating range if solution does not lie in search range

    tolerance = 1e-6 #Accuracy tolerance
    while (maxz - minz) > tolerance: #Repeat until redshift range is bigger than tolerance
        midz = (maxz + minz)/2 #calculates midpoint of redshift range
        midt = age(midz, H0, Om0, Ode0) #calculates midpoint of time range from midz and age() function
        if midt < t: #If the midpoint of the time range is less than target time, we have to 'roof' the red shift range here as red shift values higher than this will have a lower time value
            maxz = midz
        elif midt > t: #similarly, if midpoint of time range is greater than target time, the redshift range is 'floored' here as smaller values of redshift will return greater values of t
            minz = midz
        else:
            return midz #if midt = t, then mid z is the solution!
    return (maxz + minz)/2 #finally return midpoint of narrowed redshift range as solution


plt.ion()

frames = 120
num_galaxies = 60
t_range = np.linspace(0.01, 26, frames) #Simulates from 0.01 Gyr to 26 Gyr after the big bang
sf_init = 0.01

#Establishing the simulated models through a list
model1 = [70, 0.315, 0.685] #Hubble constant, Relative matter density, Relative dark energy density 
age1 = age(0, model1[0], model1[1], model1[2]) #Calculates present age of the universe (z = 0)
model2 = [70, 1, 0]
age2 = age(0, model2[0], model2[1], model2[2]) 
    
def calculate_sf(cosmo, t_value, past): #uses the z_at_time function to estimate scale factor (size) of universe at a time t_value
    z_min_search = -10000 #Initializing redshift search range
    z_max_search = 10000
    if past == False: #If simulating future, redshift search range is -1 < z < 0
        z_min_search = -0.99999
        z_max_search = 0.0
    else: #If simulating past, redshift search range is 0 < z < 10000
        z_min_search = 0.0
    
    z_solution = z_at_time(t_value, z_min_search, z_max_search, cosmo[0], cosmo[1], cosmo[2])
        
    a = 1.0 / (z_solution + 1.0) #Calculates size of universe using redshift-scale factor relation
    
    return a 

#Setting up of axis and plotting environment
fig, ax = plt.subplots(figsize=(8.5, 8))
fig.set_facecolor('black')
ax.set_facecolor('black')
ax.set_xlabel(r'X ($10^9$ Mpc)')
ax.spines['bottom'].set_color('white')
ax.spines['left'].set_color('white')
ax.set_ylabel(r'Y ($10^9$ Mpc)')
ax.set_xlim(-50, 50)
ax.set_ylim(-50, 50)
ax.set_aspect('equal', adjustable='box')
ax.tick_params(axis='x', colors='white')
ax.tick_params(axis='y', colors='white')
ax.xaxis.label.set_color('white')
ax.yaxis.label.set_color('white')
ax.grid(True, color='#444444', linestyle=':', linewidth=0.5)
age_text = ax.text(-48, 45, 'Age of Universe = 0 Gyr', color = 'white', ha='left', fontsize=12)
sf_text = ax.text(-48, 41, f'Scale Factor a(t) = {sf_init:.2f}', color = 'white', ha='left', fontsize=12)
age_text2 = ax.text(48, 45, 'Age of Universe = 0 Gyr', color = 'white', ha='right', fontsize=12)
sf_text2 = ax.text(48, 41, f'Scale Factor a(t) = {sf_init:.2f}', color = 'white', ha='right', fontsize=12)


def generate_coords(model, n): #Uses polar coordinates to generate random coordinates of galaxies for each model
    r = 45.6 * np.sqrt(np.random.uniform(0, 1, n)) #r (radius) value of random coordinate 
    if model == 1: #random theta (angle) value range determined by which model the random coordinate is generated for
        theta = np.random.uniform(np.pi*1/2, np.pi*3/2, n)
    else:
        theta = np.random.uniform(-np.pi*1/2, np.pi*1/2, n)

    #Conversion to cartesian coordinate
    x = r*np.cos(theta) 
    y = r*np.sin(theta)
    return np.vstack([x, y]).T #returns array containing 2 columns of coordinates, x and y

#Executes generation of random coordinates for both models
comoving_coords1 = generate_coords(1, int(num_galaxies/2))
comoving_coords2 = generate_coords(2, int(num_galaxies/2))

#Generates random sizes for each galaxy/point in each model
sizes1 = np.random.uniform(1,50,int(num_galaxies/2))
sizes2 = np.random.uniform(1,50,int(num_galaxies/2))

#Produces scatter object containing coordinates of all galaxies in each model
scatter1 = ax.scatter(comoving_coords1[:,0], comoving_coords1[:,1], color='yellow', marker='o', s=sizes1)
scatter2 = ax.scatter(comoving_coords2[:,0], comoving_coords2[:,1], color='yellow', marker='o', s=sizes2)

#Defines and produces arcs representing boundaries of both universes
boundary1 = Arc((0,0), width=float(91.2*sf_init), height=float(91.2*sf_init), angle=0, theta1=90, theta2=270, ec = 'red', fc = 'none', lw = 2)
boundary2 = Arc((0,0), width=float(91.2*sf_init), height=float(91.2*sf_init), angle=0, theta1=270, theta2=450, ec = 'red', fc = 'none', lw = 2)
ax.add_patch(boundary1)
ax.add_patch(boundary2)

def plot(sf, sf2): #Function runs in each frame to re-scale scatter object of galactic points and boundary objects using scale factor
    boundary1.set_width(float(91.2*sf))
    boundary1.set_height(float(91.2*sf))
    boundary2.set_width(float(91.2*sf2))
    boundary2.set_height(float(91.2*sf2))
    physical_coords1 = comoving_coords1 * sf
    physical_coords2 = comoving_coords2 * sf2
    scatter1.set_offsets(physical_coords1)
    scatter2.set_offsets(physical_coords2)
        
try:
    while True: #Indefinite while loop iterates the simulation froever
        p1 = True #Establishes past simulation of both universes initially
        p2 = True
        for i in range(frames):
            if t_range[i] > age1 and p1: #If model one begins simulating future, formatting is modifyed and a dashed arc is added to represent present size of universe (a=1)
                p1 = False
                age_text.set_color('cyan')
                sf_text.set_color('cyan')
                pboundary1 = Arc((0,0), width=float(91.2), height=float(91.2), angle=0, theta1=90, theta2=270, ec = 'cyan', fc = 'none', linestyle = '--', lw = 2)
                ax.add_patch(pboundary1)
            if t_range[i] > age2 and p2: #When model two begins simulating future
                p2 = False
                age_text2.set_color('cyan')
                sf_text2.set_color('cyan')
                pboundary2 = Arc((0,0), width=float(91.2), height=float(91.2), angle=0, theta1=270, theta2=450, ec = 'cyan', fc = 'none', linestyle = '--', lw = 2)
                ax.add_patch(pboundary2)
            if t_range[i] > np.min([age1, age2]): #If either model begins simulating the future, the axis is given larger limits
                new_axes_limit = 120
                ax.set_xlim(-new_axes_limit, new_axes_limit)
                ax.set_ylim(-new_axes_limit, new_axes_limit)
                age_text.set_position((-new_axes_limit+3, new_axes_limit-6))
                sf_text.set_position((-new_axes_limit+3, new_axes_limit-14))
                age_text2.set_position((new_axes_limit-3, new_axes_limit-6))
                sf_text2.set_position((new_axes_limit-3, new_axes_limit-14))

            #Scale factor for this frame calculated
            scale_factor = calculate_sf(model1, t_range[i], p1)
            scale_factor2 = calculate_sf(model2, t_range[i], p2)

            #Objects in axis scaled with scale factor
            plot(scale_factor,scale_factor2)
            
            #Labels updated in each frame
            age_text.set_text(f'Age of Universe = {t_range[i]:.2f} Gyr') 
            age_text2.set_text(f'Age of Universe = {t_range[i]:.2f} Gyr')
            ax.set_title(f'Standard Model (O_m0 = {model1[1]:.2f}, O_de0 = {model1[2]:.2f})    Einstein-de-sitter Model (O_m0 = {model2[1]:.2f}, O_de0 = {model2[2]:.2f})', color = 'white')
            sf_text.set_text(f'Scale Factor a(t) = {scale_factor:.2f}')
            sf_text2.set_text(f'Scale Factor a(t) = {scale_factor2:.2f}')
            fig.canvas.draw()
            plt.pause(0.05) #50ms delay between frames
        plt.pause(5) #5s delay between each full run through of simulation

        #Resetting axis and labels between each full run through of simulation
        ax.set_xlim(-50, 50)
        ax.set_ylim(-50, 50)
        pboundary1.remove()
        pboundary2.remove()
        age_text.set_color('white')
        sf_text.set_color('white')
        age_text2.set_color('white')
        sf_text2.set_color('white')
        age_text.set_position((-48, 45))
        sf_text.set_position((-48, 41))
        age_text2.set_position((48, 45))
        sf_text2.set_position((48, 41))
except KeyboardInterrupt: #End indefinite loop through keyboard interrupt
    print('Terminated')
            
plt.ioff()  
plt.show()  
