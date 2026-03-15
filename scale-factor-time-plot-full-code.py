#importing libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

#Identical functions used in simulation program
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

a_range = np.linspace(0.01, 3.0, 200) #Plot range is 0 < a(t) ≤ 3 with 200 plotting points in between

# 1. Standard Model: ~68.5% Dark Energy (DE)
cosmo_standard = [70, 0.315, 0.685]
# 2. High Dark Energy Model: 85% DE
cosmo_high_DE = [70, 0.15, 0.85]
# 3. No Dark Energy Model: 0% DE
cosmo_no_DE = [70, 1, 0]

plt.figure(figsize=(10, 6)) #Defining plot dimensions

def calculate_time(cosmo, a_values):
    times = [] #Initializes list containing x coordinates of all plots (time at which a universe has a scale factor, a)
    for a in a_values: #Iterates through all scale factor plot points 
        z = (1.0 / a) - 1.0 #Calculates redshift from scale factor using redshift-scale-factor relations
        time_at_z = age(z, cosmo[0], cosmo[1], cosmo[2]) #Evaluates time at cosmological redshift, z, for corresponding scale factor
        times.append(time_at_z) #Appends result to a list of time values
    return times

t_standard = calculate_time(cosmo_standard, a_range) #Runs function to return list of time (x) values corresponding to when the scale factors in a_range (y) are reached
plt.plot(t_standard, a_range, label=r'Standard $\Omega_{\Lambda}=0.685$', color='blue') #Plots scale factor against calculated time values

t_high = calculate_time(cosmo_high_DE, a_range)
plt.plot(t_high, a_range, label=r'High Dark Energy $\Omega_{\Lambda}=0.85$', color='red')

t_no_DE = calculate_time(cosmo_no_DE, a_range)
plt.plot(t_no_DE, a_range, label=r'Einstein-de-Sitter (No DE) $\Omega_{\Lambda}=0.0$', color='black', linestyle='--')

#Add labels and formatting
plt.title(r'Effect of Dark Energy ($\Omega_{\Lambda}$) on Cosmic Expansion ($a(t)$)', fontsize=14)
plt.xlabel('Cosmic Time (Gyr)', fontsize=12)
plt.ylabel('Scale Factor $a(t)$', fontsize=12)

# Mark the transition to acceleration for the Standard model
plt.scatter([t_standard[np.argmin(np.abs(a_range - 0.58))] ], [0.58], color='blue', marker='o', s=50, zorder=5)

#Further formatting
plt.axhline(y=1.0, color='gray', linestyle=':', linewidth=1)
plt.legend(loc='lower right')
plt.grid(True, linestyle=':', alpha=0.6)
plt.show()
