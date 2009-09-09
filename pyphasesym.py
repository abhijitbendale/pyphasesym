###############################################################################
# This module calculates the phase symmetry of points in an image.	      #
# This is a contrast invariant measure of symmetry.  This function can be     #
# used as a line and blob detector.  The greyscale 'polarity' of the lines    #
# that you want to find can be specified. The convolutions are done with FFT  #
# usage <input_filename> <output_filename>                                    #
# This module is part of pyphasesym package                                   #
###############################################################################

"""This module computes phase symmetry of points in an image"""

import optparse
import Image
import ImageOps

import numpy as np

MODULE_DTYPE = "float64"

#Default values. These are used if user does not specify these values
DEFAULT_NSCALE = 5.
DEFAULT_ORIENT = 6.
DEFAULT_MWAVELENGTH = 3.
DEFAULT_MULT = 2.1
DEFAULT_SIGMAONF = 0.55
DEFAULT_DTHETASEGMA = 1.2
DEFAULT_NSTD = 2.
DEFAULT_POLARITY = 0

class phasesym_error(Exception): 
    """ phasesym package error exception class"""
    pass

class out_of_range_error(phasesym_error): 
    """ out of range error exception"""
    pass

class array_dim_mismatch_error(phasesym_error): 
    """array mismatch error exception"""
    pass

#-------------------------------------------------------------------------------
def filter_intialization(rows, cols):

    """Filter initialization components 

    Input: no of rows and cols in the image
    Output: np arrays with filter components sintheta, costheta and radius"""
    
    # Set up X and Y matrices with ranges normalised to +/- 0.5
    # The following code adjusts things appropriately for odd and even values
    # of rows and columns.
#    assert rows > 1
#    assert cols > 1
    if not rows > 1:
        raise out_of_range_error, "cannot have image with rows < 1"
    if not cols > 1:
        raise out_of_range_error, "cannot have image with cols < 1"

    if np.mod(cols, 2):
        xr_arr = np.arange(-(cols - 1)/2, (cols - 1)/2 + 1, 
                            dtype = "float32" )/ (cols - 1)
    else:
        xr_arr = np.arange(-cols/2, (cols/2 - 1) + 1, \
                           dtype = "float32")/ cols
        
    if np.mod(rows, 2):
        yr_arr = np.arange(-(rows - 1)/2, (rows - 1)/2 + 1, \
                             dtype="float32")/ (rows-1)
    else:
        yr_arr = np.arange(-rows/2, (rows/2 - 1) +1, \
                             dtype = "float32")/ rows

    x_arr, y_arr = np.meshgrid(xr_arr, yr_arr)
    
    # Matrix values contain *normalised* radius from centre.
    radius = np.sqrt(x_arr ** 2 + y_arr ** 2)        

    # Matrix values contain polar angle.
    theta = np.arctan2(-y_arr, x_arr)

    # Quadrant shift radius and theta so that filters	       
    # are constructed with 0 frequency at the corners.       
    # Get rid of the 0 radius value at the 0		       
    # frequency point (now at top-left corner)	       
    # so that taking the log of the radius will 	       
    # not cause trouble.				       

    radius = np.fft.ifftshift(radius)
    theta = np.fft.ifftshift(theta)
    radius[0][0] = 1

    sintheta = np.sin(theta)
    costheta = np.cos(theta)

    imshape = rows, cols

    assert sintheta.shape == imshape
    assert costheta.shape == imshape
    assert radius.shape == imshape
    
    return sintheta, costheta, radius

#-------------------------------------------------------------------------------
def get_low_pass_filter(rows, cols, cutoff, n_order):

    """Compute low pass filter with given cutoff
    
    Input: imagesize, cuttoff frequency, order of filter
    Output: lowpass filter with given parameters. Implemented as np array"""

    assert 0 < cutoff <= 0.5
    assert n_order > 1
    assert rows > 1
    assert cols > 1

    if np.mod(cols, 2):
        xr_arr = np.arange(-(cols - 1)/2, (cols - 1)/2 + 1, dtype="float32")/ \
            (cols - 1)
    else:
        xr_arr = np.arange(-cols/2, (cols/2 - 1) + 1, dtype="float32")/cols

    if np.mod(rows, 2):
        yr_arr = np.arange(-(rows - 1)/2, (rows - 1)/2 + 1, dtype="float32")/ \
            (rows-1)
    else:
        yr_arr = np.arange(-rows/2, (rows/2 - 1) + 1, dtype="float32")/rows

    x_arr, y_arr = np.meshgrid(xr_arr, yr_arr)
    radius = np.sqrt(x_arr ** 2 + y_arr ** 2)
    f_comp = np.fft.ifftshift( 1./ (1. + (radius/ cutoff) ** (2*n_order) ))

    return f_comp

#-------------------------------------------------------------------------------
def get_gabor(radius, lp_filter, nscale, min_wave_length, mult, sigma_on_f):

    """ Get gabors with logarithmic transfer fucntions """
    
    log_gabor = []

    for s_count in range(int(nscale)):
        wavelength =  min_wave_length * (mult ** s_count)
        fo_comp = 1./wavelength             # Filter Frequence

        # Apply low pass filter
        log_gabor += [(np.exp((-(np.log(radius/fo_comp)) ** 2)/\
                                 (2 * np.log(sigma_on_f) ** 2))) * lp_filter]

        # Set the value at the 0 frequency point of the filter
        # back to zero (undo the radius fudge)
        log_gabor[s_count][0][0] = 0

    assert len(log_gabor) == nscale
    return log_gabor

#-------------------------------------------------------------------------------
def get_spread(sintheta, costheta, norient, d_theta_sigma):

    """ Compute angular components of the filter

    Input: sintheta, costheta (meshgrid of sin and cos values of orientation
    norient = no of orientations
    d_theta_sigma = ratio of angular interval between filter orientations
    and the std deviation of the angular Gaussian function used to construct 
    filters in the frequency plane
    Output: spread = angular component of the filter"""

    spread = []
    
    # Calculate the standard deviation of the angular Gaussian function
    # used to construct filters in the frequency plane.     
    assert norient >= 1
    assert d_theta_sigma > 0

    theta_sigma = np.pi/norient/d_theta_sigma

    # For each orientation
    for o_count in range(int(norient)):
        angl = (o_count * np.pi)/ norient   # Filter angle

        # For each point in the filter matrix calculate angular distance from
        # the specified filter orientation.  To overcome the angular wrap-around
        # problem sine difference & cosine difference values are first computed
        # and then the atan2 function is used to determine angular distance.
        
        # Difference in sine
        ds_arr = sintheta * np.cos(angl) - costheta * np.sin(angl)  
        # Difference in cosine
        dc_arr = costheta * np.cos(angl) + sintheta * np.sin(angl)  
        # Abs angular distance
        dtheta = abs(np.arctan2(ds_arr, dc_arr))                         
        # Angular filter component
        spread += [np.exp((-dtheta ** 2)/(2 * theta_sigma ** 2))]      

    assert len(spread) == norient
    return spread

#-------------------------------------------------------------------------------
def get_orientation_energy(orientation_energy, eo_conv, 
                           o_count, polarity, nscale):

    """Calculate the phase symmetry measure based on polartiy

    if polarity == 0, you get black and white spots
    if polarity == 1, you get white spots only
    if polarity ==-1 you get black spots only"""

    #look for 'white' and 'black' spots
    if polarity == 0: 
        for s_count in range(int(nscale)):
            orientation_energy = orientation_energy + \
                abs(np.real(eo_conv[s_count][o_count])) - \
                abs(np.imag(eo_conv[s_count][o_count]))

    #Just look for 'white' spots
    elif polarity == 1:
        for s_count in range(int(nscale)):
            orientation_energy = orientation_energy + \
                np.real(eo_conv[s_count][o_count]) - \
                abs(np.imag(eo_conv[s_count][o_count]))

    #Just look for 'black' spots
    elif polarity == -1:
        for s_count in range(int(nscale)):
            orientation_energy = orientation_energy - \
                np.real(eo_conv[s_count][o_count]) - \
                abs(np.imag(eo_conv[s_count][o_count]))

    return orientation_energy

#-------------------------------------------------------------------------------
def get_phasesym(rows, cols, imfft, log_gabor, 
                spread, nscale, norient, k, polarity):

    """ The main loop of phasesym 

    Input: rows, cols, imfft, log_gabor,spread, nscale, norient, k, polarity
    Output: phaseSym, orientation as numpy arrays"""

    # Array initialization
    eo_conv = np.ndarray((nscale, norient, rows, cols), dtype="complex_")
    ifft_filter_array = np.ndarray((nscale, rows, cols), dtype="float64")
    zero = np.zeros((rows, cols), dtype="float32")
    total_sum_amp = np.zeros((rows, cols), dtype="float32")
    total_energy = np.zeros((rows, cols), dtype="float32")
    orientation = np.zeros((rows, cols), dtype="float32")
    epsilon = 0.0001
    imshape = rows, cols
    
    assert norient > 0

    for o_count in range(int(norient)):  #for each orientation
        # later we combine phaswe congruency results over all orientations

        # Array initialization
        s_amp_this_orient = np.zeros((rows, cols), dtype="float32")
        orientation_energy = np.zeros((rows, cols), dtype="float32")


        for s_count in range(int(nscale)):  # for each scale
            
            #Multiply radial and angular components to get filter
            filter_comp = log_gabor[s_count] * spread[o_count]  
            ifft_filter_array[s_count] = np.real(np.fft.ifft2(filter_comp)) * \
                np.sqrt(rows * cols)
            
            # Convolve image with even and odd filters returning the result in 
            # eo_conv
            # convolving image with quadrature pair of filters
            # quadrature filter is the one where real part of the filter is 
            # related to its imaginary part via Hilbert transform along a 
            # particular axis through origin (eg. Gabor filters)
            eo_conv[s_count][o_count] = np.fft.ifft2(imfft * filter_comp)
            amp = abs(eo_conv[s_count][o_count]) # Amplitude response
            s_amp_this_orient = s_amp_this_orient + amp
            
            # Record mean squared filter value at smallest
            # scale. This is used for noise estimation
            if s_count == 0:
                em_n = sum(sum(filter_comp ** 2)) 
       
        # Now Calulate phase symmetry measure
        orientation_energy = get_orientation_energy(orientation_energy, 
                                                    eo_conv, o_count, 
                                                    polarity, nscale)


        assert orientation_energy.shape == imshape
        # Noise Compensation
        # We estimate the noise power from the energy squared response at the
        # smallest scale.  If the noise is Gaussian the energy squared will
        # have a Chi-squared 2DOF pdf.  We calculate hte median energy squared
        # response as this is a robust statistic.  From this we estimate the
        # mean.  The estimate of noise power is obtained by dividing the mean
        # squared energy value by the mean squared filter value

        median_e2n = np.median((abs(eo_conv[0][o_count]) ** 2).ravel())
        mean_e2n = -median_e2n/np.log(0.5)

        noise_power = mean_e2n/em_n           # Estimate noise power

        # Now estimate the total energy^2 due to noise
        # Estimate for sum(An^2) + sum(Ai.*Aj.*(cphi.*cphj + sphi.*sphj))

        est_sum_an_2 = np.zeros((rows, cols), dtype="float32")
        for s_count in range(int(nscale)):
            est_sum_an_2 = est_sum_an_2 + ifft_filter_array[s_count] ** 2

        est_sum_ai_aj = np.zeros((rows, cols), dtype="float32")
        for si_count in range(int(nscale - 1)):
            for sj_count in range(si_count+1, int(nscale)):
                est_sum_ai_aj = est_sum_ai_aj + \
                    ifft_filter_array[si_count] * \
                    ifft_filter_array[sj_count]
        
        est_noise_energy_2 = 2 * noise_power * sum(sum(est_sum_an_2)) \
            + 4 * noise_power * sum(sum(est_sum_ai_aj))

        tau = np.sqrt(est_noise_energy_2/2)        # Rayleigh parameter
        # Expected value of noise energy
        est_noise_energy = tau * np.sqrt(np.pi/2) 
        est_noise_energy_sigma = np.sqrt((2 - np.pi/2) * tau**2)

        # Noise threshold
        noise_thresh =  est_noise_energy + k * est_noise_energy_sigma  

        # The estimated noise effect calculated above is only valid for the PC_1
        # measure.  The PC_2 measure does not lend itself readily to the same
        # analysis.  However empirically it seems that the noise effect is
        # overestimated roughly by a factor of 1.7 for the filter parameters
        # used here.
        noise_thresh = noise_thresh/1.7

        # Apply noise threshold
        orientation_energy = np.maximum(orientation_energy - noise_thresh, zero)

        # Update accumulator matrix for sumAn and total_energy
        total_sum_amp = total_sum_amp + s_amp_this_orient
        total_energy = total_energy + orientation_energy

        # Update orientation matrix by finding image points where the energy in
        # this orientation is greater than in any previous orientation (the
        # change matrix) and then replacing these elements in the orientation
        # matrix with the current orientation number.

        if o_count == 0:
            max_energy = orientation_energy
        else:
            change = orientation_energy > max_energy
            orientation = o_count * change + orientation * \
                np.logical_not(change)
            max_energy = np.maximum(max_energy, orientation_energy)

    # Normalize total_energy by the total_sum_amp to obtain phase symmetry
    # epsilon is used to avoid division by 000
    phasesym_arr = total_energy / (total_sum_amp + epsilon)
    assert phasesym_arr.shape == imshape

    # Convert orientation matrix values to degrees
    orientation = orientation * (180/norient)
    assert orientation.shape == imshape
        
    return phasesym_arr, orientation

#-------------------------------------------------------------------------------
def phasesym(input_array, nscale, norient, min_wave_length, mult, sigma_on_f,
             d_theta_sigma, k, polarity):

    """ Modular interface for various operations in phasesym

    Input: image as np array and other arguments for phasesym. Check
    Readme for further details about arguments. If no arguments are
    provided, the code uses default arguments. However, image as numpy array
    must be provided
    Output: phaseSym and orientation of image as numpy arrays"""
    
    rows, cols = input_array.shape

    imfft = np.fft.fft2(input_array)

    assert input_array.shape == imfft.shape

    # Filter initializations
    sintheta, costheta, radius = filter_intialization(rows, cols)

    # Filters are constructed in terms of two components.
    # 1) The radial component, which controls the frequency band that the filter
    #    responds to
    # 2) The angular component, which controls the orientation that the filter
    #    responds to.
    # The two components are multiplied together to construct the overall filter

    # Construct filter radial components
    
    # First construct a low-pass filter that is as large as possible, yet falls
    # away to zero at the boundaries.  All log Gabor filters are multiplied by
    # this to ensure no extra frequencies at the 'corners' of the FFT are
    # incorporated as this seems to upset the normalisation process when
    # calculating phase congrunecy.
    
    # Construct Low pass filter
    lp_filter = get_low_pass_filter(rows, cols, 0.4, 10.)
    
    assert input_array.shape == lp_filter.shape

    # Radial Component
    log_gabor = get_gabor(radius, 
                          lp_filter,
                          nscale,
                          min_wave_length,
                          mult,
                          sigma_on_f)

    # Construct the angular filter components
    spread = get_spread(sintheta,
                        costheta,
                        norient,
                        d_theta_sigma)

    # Get phase symmetry and orientation of image
    phasesym_arr, orientation = get_phasesym(rows, 
                                             cols, 
                                             imfft,
                                             log_gabor,
                                             spread,
                                             nscale,
                                             norient,
                                             k,
                                             polarity)

    return phasesym_arr, orientation

#-------------------------------------------------------------------------------
def phasesym_from_array(input_array, nscale, norient, min_wave_length, mult,
                       sigma_on_f, d_theta_sigma, k, polarity):

    """ Calculate phasesym from image as numpy array

    Inputs: image as an np array, and other parameters for phasesym. 
    input_array = image as a numpy array
    nscale = Number of wavelet scales, try values 3-6		        
    norient = Number of filter orientations.			        
    min_wave_length = Wavelength of smallest scale filter.		    
    mult = Scaling factor between successive filters.		    
    sigma_on_f = ratio of std deviation of gaussians for logGabors
    d_theta_sigma = Ratio of angular interval between filter orientations  
    and the standard deviation of the angular Gaussian	    
    function used to construct filters in the		    
    freq. plane.					    
    k = No of standard deviations of the noise energy 
    polarity = (0, 1, -1) for bright and dark points, bright points
    or dark points respectively

    Outputs: numpy arrays 
    phaseSym = phase symmetry of image
    Orientation = Orientation Image
    """
    
    # Call to phasesym
    phasesym_arr, orientation = phasesym(input_array, 
                                         nscale, 
                                         norient, 
                                         min_wave_length, 
                                         mult, 
                                         sigma_on_f, 
                                         d_theta_sigma, 
                                         k, 
                                         polarity)

    assert input_array.shape == phasesym_arr.shape
    assert input_array.shape == orientation.shape

    return phasesym_arr, orientation

#-------------------------------------------------------------------------------
def phasesym_from_filename(
    input_filename,
    output_filename,
    nscale,
    norient,
    min_wave_length,
    mult,
    sigma_on_f,
    d_theta_sigma,
    k,
    polarity
    # Here we can also have check for output fileformat
    ):
    
    """Basic input output file handling function.
    
    Input: input and output filenames and phasesym arguments
    This function also invokes call to phasesym_fromArray function
    This function also saves the output from phasesym calculation
    in user-specified format and path. For more details on output format
    check README"""    

    # Read input Image
    img = Image.open(input_filename)
    img = ImageOps.grayscale(img)
    imarr = np.asarray(img)    

    # Call to phasesym
    phasesym_arr, orientation = phasesym_from_array(imarr, 
                                                    nscale, 
                                                    norient, 
                                                    min_wave_length,
                                                    mult, 
                                                    sigma_on_f, 
                                                    d_theta_sigma, 
                                                    k, 
                                                    polarity)
    
    assert imarr.shape == phasesym_arr.shape
    assert imarr.shape == orientation.shape

    print phasesym_arr
    print orientation
    # pkl
    import cPickle
    
    data = {'phaseSym': phasesym_arr,
            'orientation': orientation,
            }
    cPickle.dump(data, open(output_filename, "w+"), protocol=2)

    
#-------------------------------------------------------------------------------
def main():

    """ Main() function and optparsing """

    usage = "usage: %prog [options] <input_filename> <output_filename>"

    parser = optparse.OptionParser(usage=usage)

    parser.add_option("--nscale", "-s",
                      default=DEFAULT_NSCALE,
                      type="float", 
                      metavar="FLOAT",
                      help="[default=%default]")
    
    parser.add_option("--norient", "-o",
                      default=DEFAULT_ORIENT,
                      type="float", 
                      metavar="FLOAT",
                      help="[default=%default]")

    parser.add_option("--min_wave_length", "-w",
                      default=DEFAULT_MWAVELENGTH,
                      type="float", 
                      metavar="FLOAT",
                      help="[default=%default]")

    parser.add_option("--mult", "-m",
                      default=DEFAULT_MULT,
                      type="float", 
                      metavar="FLOAT",
                      help="[default=%default]")

    parser.add_option("--sigma_on_f", "-g",
                      default=DEFAULT_SIGMAONF,
                      type="float", 
                      metavar="FLOAT",
                      help="[default=%default]")

    parser.add_option("--d_theta_sigma", "-d",
                      default=DEFAULT_DTHETASEGMA,
                      type="float", 
                      metavar="FLOAT",
                      help="[default=%default]")

    parser.add_option("--nstdeviations", "-k",
                      default=DEFAULT_NSTD,
                      type="float", 
                      metavar="FLOAT",
                      help="[default=%default]")

    parser.add_option("--polarity", "-p",
                      default=DEFAULT_POLARITY,
                      type="int", 
                      metavar="INT",
                      help="[default=%default]")

    opts, args = parser.parse_args()


    if len(args) < 1:
        print "ERROR: Supply Image"
        parser.print_help()
    else:
        input_filename = args[0]
        output_filename = args[1]
        
        phasesym_from_filename(input_filename,
                              output_filename,
                              # -- 
                              nscale=opts.nscale,
                              norient=opts.norient,
                              min_wave_length=opts.min_wave_length,
                              mult=opts.mult,
                              sigma_on_f=opts.sigma_on_f,
                              d_theta_sigma=opts.d_theta_sigma,
                              k=opts.nstdeviations,
                              polarity=opts.polarity
                              )


#-------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
