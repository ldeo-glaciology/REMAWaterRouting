# destripe an image with oblique (i.e. angled stripes)

import numpy as np
import matplotlib.pyplot as plt
import cv2
from Oblique_destriping_functions import striping_angle, grad_image,\
super_Gauss_filter, refinement_destriping

striped_image = np.load('LS8_patch.npy')
Ly_i,Lx_i = np.shape(striped_image)


theta = striping_angle(striped_image)

if np.amax(grad_image(striped_image/np.mean(striped_image)))<0.006:
    destriped_image = super_Gauss_filter(striped_image,0.49*Lx_i,110,theta)
    
else:
    N_pad = 300
    padded_patch = np.zeros((Ly_i+2*N_pad,Lx_i+2*N_pad))
    padded_patch[N_pad:N_pad+Ly_i,N_pad:N_pad+Lx_i] = striped_image
    
    Ly_p = Ly_i + 2*N_pad
    Lx_p = Lx_i + 2*N_pad
    
    M_forward = cv2.getRotationMatrix2D((Lx_p/2,Ly_p/2),theta,1)
    array_rot = cv2.warpAffine(padded_patch,M_forward,(Lx_p,Ly_p))
    
    padded_array_rot = np.copy(array_rot)
    padded_array_rot[array_rot==0] = np.mean(striped_image[striped_image>0])
    
    # use standard values for filter parameters
    destriped_array_rot = super_Gauss_filter(padded_array_rot,Lx_p*0.485,100,0)
    
    # refine guess of stripe contribution to image
    destriped_array_ref = refinement_destriping(array_rot,destriped_array_rot)
    
    M_backward = cv2.getRotationMatrix2D((Lx_p/2,Ly_p/2),-theta,1)
    destriped_array_unrot = cv2.warpAffine(destriped_array_ref ,M_backward,\
                                           (Lx_p,Ly_p))
    destriped_image = destriped_array_unrot[N_pad:N_pad+Ly_i,N_pad:N_pad+Lx_i]
    
    
    
    #%%
plt.rcParams['figure.figsize'] = [8, 4]    
plt.subplot(1,2,1)
plt.imshow(striped_image)
plt.colorbar()
plt.clim(24300,25000)

plt.subplot(1,2,2)
plt.imshow(destriped_image)
plt.colorbar()
plt.clim(24300,25000)


