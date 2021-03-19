import numpy as np
import cv2

# utilities

# extract angle of oblique stripes
def stripe_angle(array):
    FT = np.fft.fftshift(np.fft.fft2(array))
    Ly,Lx = np.shape(array)
    LO2 = np.fft.fftshift(np.abs(FT[:,0])) # take lineout from first column
    #peaks, _ = find_peaks(LO2, height=15.5,distance=15) # find location of peak from strip
    peaks2 = np.argmax(LO2[int(np.round(Lx/12)):int(np.round(Lx/2))])
    peaks2 = peaks2 + int(np.round(Lx/12))
    # calculate angle of stripe
    theta_in  = np.rad2deg(np.arctan((2*(peaks2+1)/Ly)))# convert to degrees for simplicity
    return theta_in

def striping_angle(array):
    FT = np.fft.fftshift(np.fft.fft2((array-np.mean(array))))
    Ly,Lx = np.shape(FT)
    
    # scan quadrant of Fourier space for sidband position
    SB = []
    for ii in range(0,int(Lx/2)):
        LO = np.abs(FT[:int(Ly/2),ii]) # lineout
        # assume sideband is peak of lineout
        SB.append(np.argmax(LO))
    
    # initial linear fit to sideband locations to help remove outliers
    test_vec = np.asarray(SB)
    xx = np.arange(0,int(Lx/2),1)
    coefs = np.polyfit(xx,test_vec,1)
    
    # exclude outliers based on first fit
    adj_vec = np.asarray(SB - (coefs[0]*xx+coefs[1]))
    
    coefs_adj = np.polyfit(xx[np.abs(adj_vec-np.mean(adj_vec))<np.std(adj_vec)]\
                          ,test_vec[np.abs(adj_vec-np.mean(adj_vec))<np.std(adj_vec)],1)
    
    # actual sideband angle related to gradiet of fit to sideband locations
    theta=  np.degrees( np.arctan(coefs_adj[0]))
    
    return theta


def grad_image(array):
    return (cv2.Laplacian(np.asarray(array,dtype = np.float32),cv2.CV_32FC1))

# functions used for destriping

# Rotated super-Gaussian for the filter
def SGaussian_rotated(x_mesh,y_mesh,x0,y0,wx,wy,px,py,theta):
    theta = np.deg2rad(theta)
  
    x_mesh_rot = x_mesh*np.cos(theta)-y_mesh*np.sin(theta)
    y_mesh_rot = x_mesh*np.sin(theta)+y_mesh*np.cos(theta)
    
    y_mesh_rot = np.fft.fftshift(y_mesh_rot,1)   
    
    SG = np.exp(-abs((x_mesh_rot-x0)/wx)**px-abs((y_mesh_rot-y0)/wy)**py)
   
    SG = 1-SG
    return SG


def super_Gauss_filter(input_img,sg_width,sg_pow,theta):
    mask = np.copy(input_img) # copy avoids writing issues

    Ly,Lx = np.shape(input_img)
    x = np.linspace(-Lx/2-1, Lx/2,Lx)
    y = np.linspace(-Ly/2-1, Ly/2,Ly)

    x_mesh, y_mesh = np.meshgrid(x, y)


    # generate filter and apply to Fourier transform
    #recall SG(x_mesh,y_mesh,offset in x, offset in y, width across Fourier domain band,
    # width along Fourier domain band, super-Gaussian power across band, angle in Fourier domain)
    SG = SGaussian_rotated(x_mesh,y_mesh,0,0,Lx/200,sg_width,2,sg_pow,90-theta) # 90 degree angle for verical stripes

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
    #output_vec[(np.abs(output_vec)-np.mean(output_vec))>2*np.std(output_vec)]=0 

    w1 = np.ones(Ly)
    stripe_matrix = np.outer(w1,output_vec)
    refined_destriped_image = striped_image-stripe_matrix
    
    return refined_destriped_image
