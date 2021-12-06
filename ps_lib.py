from pathlib import Path
from typing import Union

import numpy as np
import numpy.random
from PIL import Image
import matplotlib.colors as mc
__all__ = ['read_image', 'write_image', 'angular_fourier_frequencies']

#------------------------------------------------------------------------------

def read_image(path: Union[Path, str]) -> np.ndarray:
    '''
    Read a PNG or JPG image an array of linear RGB radiance values ∈ [0,1].
    '''
    return (np.float32(Image.open(path)) / 255)**2.2


def write_image(path: Union[Path, str], image: np.ndarray, alpha = 1.0, jh = 0., js = 0., jv = 0., jr = 0., jg = 0., jb = 0.) -> None:
    '''
    Write an array of linear RGB radiance values ∈ [0,1] as a PNG or JPG image.
    '''
    '''
    # Jitter the hsv image as needed
    if jh != 0. or js != 0. or jv != 0.:
        hsvimage = mc.rgb_to_hsv(image)
        
        # Perturb as specified
        hsvimage[:,:,0] += np.random.uniform(-jh, jh, hsvimage[:,:,0].shape)
        
        hsvimage[:,:,1] += np.random.uniform(-js, js, hsvimage[:,:,1].shape)
        
        hsvimage[:,:,2] += np.random.uniform(-jv, jv, hsvimage[:,:,2].shape)
        
        # Transfer back to rgb to add opacity
        image_jittered = mc.hsv_to_rgb(hsvimage.clip(0, 1))
    
    # Jitter the rgb image
    else:
        image_jittered = np.copy(image)
        # Perturb as specified
        image_jittered[:,:,0] += np.random.uniform(-jr, jr, image_jittered[:,:,0].shape)
        image_jittered[:,:,1] += np.random.uniform(-jg, jg, image_jittered[:,:,1].shape)
        image_jittered[:,:,2] += np.random.uniform(-jb, jb, image_jittered[:,:,2].shape)
        
    '''
    # Jitter the hsv image as needed
    if jh != 0. or js != 0. or jv != 0.:
        hsvimage = mc.rgb_to_hsv(image)
        
        # Perturb as specified
        hsvimage[:,:,0] *= (1 + np.random.uniform(-jh, jh))
        
        hsvimage[:,:,1] *= (1 + np.random.uniform(-js, js))
        
        hsvimage[:,:,2] *= (1 + np.random.uniform(-jv, jv))
        
        # Transfer back to rgb to add opacity
        image_jittered = mc.hsv_to_rgb(hsvimage.clip(0, 1))
    
    # Jitter the rgb image
    else:
        image_jittered = np.copy(image)
        # Perturb as specified
        image_jittered[:,:,0] *= (1 + np.random.uniform(-jr, jr))
        image_jittered[:,:,1] *= (1 + np.random.uniform(-jg, jg))
        image_jittered[:,:,2] *= (1 + np.random.uniform(-jb, jb))
    my_image = np.zeros([image_jittered.shape[0], image_jittered.shape[1], 4], dtype='uint8')
    
    # Add opacity
    my_image[:,:,0:3] = np.uint8(255 * image_jittered.clip(0, 1)**(1/2.2))
    my_image[:,:,3] = np.uint8(255 * np.clip(alpha, 0, 1))
    
    # print(my_image[:,:,3])
    Image.fromarray(np.uint8(my_image), mode = 'RGBA').save(path)

def write_image_Alpha(path: Union[Path, str], image: np.ndarray) -> None:
    '''
    Write an array of linear RGB radiance values ∈ [0,1] as a PNG or JPG image.
    '''
    
    
    # Add opacity
    my_image = np.uint8(255 * image.clip(0, 1)**(1/2.2))
    
    # print(my_image[:,:,3])
    Image.fromarray(np.uint8(my_image), mode = 'RGBA').save(path)

def angular_fourier_frequencies(height: int, width: int) -> np.ndarray:
    '''
    Return angluar freqeuncy components for each DFT component.
    '''
    col_freqs = np.fft.fftfreq(height)[:, None]
    row_freqs = np.fft.fftfreq(width)[None, :]
    col_freqs[0, :] = col_freqs[1, :]
    row_freqs[:, 0] = row_freqs[:, 1]
    return 2 * np.pi * np.sqrt(row_freqs**2 + col_freqs**2)
