# This code is used to plot the power spectrum of the image and the wavelet power spectrum

#We've set this code up so that you have to input certain things and learn some of how we analyze data
#

#Imports to be used in script
#It is generally good practice to only import what is used in the code
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# You'll need this when making the imshow be log
from matplotlib.colors import LogNorm 


# Set range for fitting the power spectrum in pixel units. This will be used later
lmin_fit = 2
lmax_fit = 128

# These next three lines adjust the font size for axis labels and title
plt.rcParams.update({'font.size': 20})  # Set the base font size
plt.rcParams['axes.labelsize'] = 20     # Set the axis label size
plt.rcParams['axes.titlesize'] = 20     # Set the title size

# We now want to read in the data which in astronomy is typically given as a fits file
# The most common way to read in a fits file in Python is with astropy.io's fits.open command
# Below is the basic set up to read in a fits file.
#Your first task is to find the path on your computer for casa_0.5-1.5keV.fits and read it in
with fits.open() as hdul:
    # we now want the actual data in the fits file which we get from the command below
    data = hdul[0].data
    
    #Try printing the data to see what it looks like  (should be an extremly large array)
    # We want to make sure that any values less than 0 are set to 0
    # Write a line or two of code to set that. Feel free to look it up online or to ask a chatbot
    
    
    
    # Next, we get the pixel scale from the FITS header (assumed to be in degrees)
    # Fits headers contain a ton of useful information try printing one out to see all the information
    pixel_scale = abs(hdul[0].header.get('CDELT1', 1)) *60*60 # Convert degrees to arcseconds
    
    # Calculate extent of object in arcminutes
    # You will want to use this when plotting the image
    height, width = data.shape #since the data is 2D it will assign the first value of the .shape command to height and the second value to width
    extent = [-width/2 * pixel_scale, width/2 * pixel_scale, 
              -height/2 * pixel_scale, height/2 * pixel_scale]
    
    
    #Let's plot the data to see what it looks like first by initializing a figure
    plt.figure(figsize=(10, 8))
    
    # since this is an image we want to use plt.imshow
    # look up the documentation for it to see what needs to be put in to view the data
    # There will also be optional keywords we might want to try 
    # I recommend plotting the data in log scale (which imshow has a command to do just that)
    plt.imshow()

    #If everything has gone correctly, we also want to include a color bar which I've provided the code for below 
    plt.colorbar(label='Intensity')
    
    
    # See all that white on the edges? we don't really care about it because those are just random pixels
    # We want to focus on the the actual object
    # Write some lines of code that restrict the x and y axis so this is focused around what we're interested in
    
    #Optional: Add title (using object name from header if available)
    
        
    #Optional: Add axis labels
  
    plt.show()



    # We now want to look at the Wavelet power spectrum as a function of wavelet scale
    # Wavelet transform
    # Define scale range for wavelet transform using powers of 2
    # The scale is the size of the wavelet in pixels
    # Vary the range to match the image size
    scales = 2**np.arange(0, np.log2(min(height,width))-1)  # Gives [1,2,4,8,...,2048]
    # Initialize list for wavelet power
    wavelet_power = []
    
 # Calculate wavelet power spectrum
 # The wavelet is a Ricker wavelet with an L2 normalization
    for scale in scales:
            # Create Mexican hat wavelet filter
            x = np.linspace(-scale*4, scale*4, int(scale*8))
            y = x[:, np.newaxis]
            wavelet = np.exp(-(x**2 + y**2)/(2*scale**2)) * (1 - (x**2 + y**2)/(2*scale**2))
                  # L2 normalize the wavelet
            wavelet = wavelet / np.sqrt(np.sum(wavelet**2))      
            # Convolve with image
            conv = np.abs(np.fft.ifft2(np.fft.fft2(data) * np.fft.fft2(wavelet, data.shape)))
            wavelet_power.append(np.mean(conv**2))


    # Plot wavelet power spectrum
    plt.figure(figsize=(10, 8))

    # Fit power law only to scales between lmin and lmax defined earlier
    mask = (scales >= lmin_fit) & (scales <= lmax_fit)
    scales = scales*abs(hdul[0].header.get('CDELT1', 1)) *60*60

    # Your task is to plot the power versus the scale and then get a best fit slope 
    # hint think about what would be easiest to then fit a line of best fit to
    # start with something like plt.plot and then explore from there

   
    # we want to find the line of best bit for the data. To do that we use np.polyfit
    # I reccomend reading documentation on polyfit to understand how it works
    # Fill in the terms for the function so that we get a best fit slope and intercept
    # This can be tricky for your first time so feel free to ask questions!
    fit = np.polyfit()
    
    # Now plot the line of best fit on the same plot 
    # You'll need to think about how to use the slope and intercept to get the line to show
    
    
   
    
    #uncomment these few lines to have the best fit slope be reported and for a grid on the plot
    
    #plt.text(0.05, 0.95, f'Best fit slope: {fit[0]:.2f}', 
     #        transform=plt.gca().transAxes, 
      #       bbox=dict(facecolor='white', alpha=0.8))
    #plt.grid(True, which='both', linestyle='--', alpha=0.5)
    
    
    #label the x and y axis
    
   
    # Something we'll want to do with the figures is save them. 
    #If you're using Spyder you can save them individually but that takes time and can have bad formatting
    #Instead I recommend using plt.savefig
    # I've included a formatting argument in the line that makes my plots look a little nicer
    #When you're satisfied with your plot, save it as a png or pdf. 
    #You'll want to give it an easily recognizable name and save it to a place you can easily access.
    
    #plt.savefig(FILENAME, ,bbox_inches='tight')
    plt.show()
  
   
  
  