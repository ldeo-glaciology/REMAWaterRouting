# Destriping functions: striped image in, destriped image out
# assumes stripes aligned vertically in image

import numpy as np

# Guan et al libraries
from keras.models import Model
from keras.layers import  Input, Conv2D, Activation, Add

# Pystripe library
import pystripe

# Pande-Chhetri libraries
import pywt

# Rogass et al. libraries
from scipy import signal


#%% Our technique (oriented super-Gaussian filter)

import cv2
# calculate approximate image gradient (for gradient test)
def grad_image(array):
    return (cv2.Laplacian(np.asarray(array),cv2.CV_32FC1))

# Rotated super-Gaussian for the filter
def SGaussian_rotated(x_mesh,y_mesh,x0,y0,wx,wy,px,py,theta):
    theta = np.deg2rad(theta)
  
    x_mesh_rot = x_mesh*np.cos(theta)-y_mesh*np.sin(theta)
    y_mesh_rot = x_mesh*np.sin(theta)+y_mesh*np.cos(theta)
    
    y_mesh_rot = np.fft.fftshift(y_mesh_rot,1)   
    
    SG = np.exp(-abs((x_mesh_rot-x0)/wx)**px-abs((y_mesh_rot-y0)/wy)**py)
   
    SG = 1-SG
    return SG


def super_Gauss_filter(input_img,sg_width,sg_pow):
    mask = np.copy(input_img) # copy avoids writing issues

    Ly,Lx = np.shape(input_img)
    x = np.linspace(-Lx/2-1, Lx/2,Lx)
    y = np.linspace(-Ly/2-1, Ly/2,Ly)

    x_mesh, y_mesh = np.meshgrid(x, y)


    # generate filter and apply to Fourier transform
    #recall SG(x_mesh,y_mesh,offset in x, offset in y, width across Fourier domain band,
    # width along Fourier domain band, super-Gaussian power across band, angle in Fourier domain)
    SG = SGaussian_rotated(x_mesh,y_mesh,0,0,Lx/200,sg_width,2,sg_pow,90) # 90 degree angle for verical stripes

    # apply filter in Fourier domain
    FT = np.fft.fft2(mask)
    FT = np.fft.fftshift(FT)
    FT_filt = FT*(SG) 
    FT_filt = np.fft.ifftshift(FT_filt)
    img_filt = np.fft.ifft2(FT_filt)
    img_filt = abs(img_filt)
    
    # adjust filtered pixels to account for "power" loss caused by filter
    # as per Parseval's theorem
    scale_fac = np.sum(mask)/np.sum(img_filt)
    img_filt = img_filt *scale_fac

    return img_filt

# Using estimate of stripes provided after filtering, refine the estimate
# and destripe with that
def refinement_destriping(striped_image,destriped_image):
    
    stripe_estimate = striped_image - destriped_image
    Ly,Lx = np.shape(stripe_estimate)
    
    output_vec = np.zeros(Lx)
    for ii in range(0,Lx):
        column = stripe_estimate[:,ii]
        column = column[striped_image[:,ii]>0]
        
        # column length non-zero take median, else do nothing to destriping image
        if np.sum(np.abs(column))>0:
            output_vec[ii] = np.median(column)
        else:
           output_vec[ii] = 0
           
    # exclude outliers
    output_vec[(np.abs(output_vec)-np.mean(output_vec))>4*np.std(output_vec)]=0 
        
    w1 = np.ones(Ly)
    stripe_matrix = np.outer(w1,output_vec)
    refined_destriped_image = striped_image-stripe_matrix
    
    return refined_destriped_image

def Lloyd_destripe(array,width,power):
    Ly,Lx = np.shape(array)
    filtered_array = super_Gauss_filter(array,width,power)
    
    #gradient test
    if np.amax(grad_image(np.asarray(array/np.mean(array),dtype= np.float32)))>0.1:
    
        # refine guess of stripe contribution to image
        destriped_array = refinement_destriping(array,filtered_array)
        return destriped_array
    else:
        return filtered_array
     

    


#%% Guan et al. (needs weight file)

# tests indicate that image needs to be min/max normalised for best performance
def minmax_norm(array,min_val,max_val):
    return (array-min_val)/(max_val-min_val)

def undo_norm(array,min_val,max_val):
    return array*(max_val-min_val)+min_val

def SNRDWNN():

    inpt = Input(shape=(None,None,4))
    x = Conv2D(filters=64, kernel_size=(3,3), strides=(1,1), padding='same' ,kernel_initializer=init,name='Conv-1')(inpt)
    x = Activation('relu')(x)
    for i in range(8):
        x = Conv2D(filters=64, kernel_size=(3,3), strides=(1,1), padding='same' ,kernel_initializer=init)(x)
        x = Activation('relu')(x)
    residual = Conv2D(filters=4, kernel_size=(3,3), strides=(1,1), padding='same' ,kernel_initializer=init, name = 'residual')(x)
    res = Add(name = 'res')([inpt,residual])
    model = Model(inputs=inpt, 
                  outputs=[res,residual],
                  name = 'DWSRN'
                  )

    return model

L2 = None
init = 'he_normal'
checkpoint_file = 'weights'
WEIGHT_PATH = './'+checkpoint_file+'/weight.hdf5'
model =  SNRDWNN()
model.load_weights(WEIGHT_PATH)


def Guan_destripe(input_array):
                    
    array = minmax_norm(input_array,np.amin(input_array),np.amax(input_array))
    # Discrete Wavelet Transform
    LLY,(LHY,HLY,HHY) = pywt.dwt2(array, 'haar')
    Y = np.stack((LLY,LHY,HLY,HHY),axis=2)
    # predict
    x_test = np.expand_dims(Y,axis=0)
    y_pred,noise = model.predict(x_test)

    coeffs_pred = y_pred[0,:,:,0],(LHY,y_pred[0,:,:,2],HHY)

    img_out = pywt.idwt2(coeffs_pred, 'haar')

    output_array = np.clip(img_out, 0, 1) # necessary?
    
    output_array = undo_norm(output_array,np.amin(input_array),np.amax(input_array))
    return output_array

#%% Rogass et al.
    # Rogass et al.

def Rogass_destripe(input_array):
    
    rad_gradx = np.gradient(input_array,axis=1) # gradient in x (horizontal) direction

    # "boxcar filter of length 3: middle column entries have value 1/3 (for normalisation)
    fil2 = np.zeros([3,3])
    fil2[1] = 1/3
    fil2 = np.transpose(fil2)

    # gradient in x
    rad_gradx_sup = signal.convolve2d(rad_gradx,fil2,boundary='fill', mode='same')

    LHS = np.median(rad_gradx_sup,axis=0) #median of rows. 0 or 1?

    w1 = np.ones(np.shape(input_array)[0]) # "a row vector of width equal to the number of rows and valued 1" 

    # integrand of equation (6) in Rogass et al.
    integrand = np.outer(w1,LHS) 
    # not sure of the point in integrating the integrand over x, so ignore here

    #gradient calculation causes one pixel shift, roll ==matlab circshift, corrected here
    stripe_guess_1 = np.roll(integrand,1,1) # equation (6) ignoring integral


    # recommended post-processing steps

    # equation (7)
    stripe_guess_2 = stripe_guess_1-np.mean(stripe_guess_1)
    stripe_guess_2 = stripe_guess_2 - np.mean(input_array)\
                    + np.mean(input_array-stripe_guess_1+np.mean(stripe_guess_1))
    
    # skip equation (8) as no stray light or long wave gradients in this problem
    
    output_array = input_array - stripe_guess_2 # destriped image
    return output_array

#%% Pande-Chhetri

def frequency_domain_filtering(HL,k):
    
    # Initialize the filter HLf
    HLf = np.zeros(HL.shape)
    
    # Determine the number of columns
    cols = HL.shape[1]
    
    for i in range(cols):
        # Extract the column x
        x = HL[:,i]
        
        # Compute the fft of x
        Xf = np.fft.fft(x)
        
        # Compute the mean of the vector
        mu_x = np.mean(x)
        
        # Compute the standard deviation of the vector
        sig = np.std(x)
        
        # Extract the sub-vector y
        y = x[np.abs(x - mu_x) < k*sig]
        
        # Compute the mean of the vector
        mu_y = np.mean(y)
        
        # Compute the fft of y
        Yf = np.fft.fft(y)
        
        # Extract the original dc component
        Forig = Yf[0]
        
        # Derive the normalized dc component
        Fnorm = Forig * (mu_x - mu_y)/mu_x
        
        # Replace the dc component of Xf
        Xf[0] = Fnorm
        
        # Compute the invese fft
        HLf[:,i] = np.real(np.fft.ifft(Xf))

    return HLf

def Pande_Chhetri_destripe(img_in,L,k):
# Compute the multiscale wavelet decomposition db4 level 4 as per paper
    coeffs = pywt.wavedec2(img_in, wavelet='db4', level=int(L))

    # Filter the stripes in the frequency domain
    coeffs_filt = np.copy(coeffs)
    for ii in range(0,int(L)):
        coeffs_filt[ii+1][1][:] = frequency_domain_filtering(coeffs[ii+1][1],k)

    img_tilde = pywt.waverec2(tuple(coeffs_filt), 'db4')
    return img_tilde
    
#%% Pystripe
# awkwardly needs horizontal stripes and an expanded dynamic range to 
# work well
def Pystripe_destripe(array, sigma_1,sigma_2,level_):
    factor = 2**12
    striped_array = np.transpose(array)
    mod_fac = np.amax(striped_array)/factor
    striped_array = striped_array/mod_fac
    destriped_array = pystripe.filter_streaks(striped_array, sigma=[sigma_1, sigma_2],\
                                              level=level_, wavelet='db2')
    #undo modifications
    destriped_array =np.transpose(destriped_array*mod_fac)
    
    return destriped_array

#%% define some standard performance metrics
    
from skimage.metrics import structural_similarity as ssim
from skimage.metrics import peak_signal_noise_ratio as psnr

def RMSE_metric(img1,img2):
    img1 = np.asarray(img1)
    img2 = np.asarray(img2)
    return(np.sqrt(np.mean( (img1 - img2) ** 2 )))
    
def SSIM_metric(img1,img2):
    img1 = np.asarray(img1)
    img2 = np.asarray(img2)
    pixel_max = np.max(np.maximum(img1, img2))
    pixel_min = np.min(np.minimum(img1, img2))
    pixel_diff = pixel_max-pixel_min

    return(ssim(img1,img2,data_range = pixel_diff))

def PSNR_metric(img1,img2):
    img1 = np.asarray(img1)
    img2 = np.asarray(img2)
    pixel_max = np.max(np.maximum(img1, img2))
    pixel_min = np.min(np.minimum(img1, img2))
    pixel_diff = pixel_max-pixel_min
    return(psnr(img1,img2,data_range = pixel_diff))