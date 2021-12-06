#!/usr/bin/env python3

import imageio
import numpy as np
from pathlib import Path
from ps_lib import *
import scipy as sp
import scipy.ndimage
import random
#------------------------------------------------------------------------------
# Define an answer generator for each question.

def answer_q1() -> None:
    # Load the image
    imagename = 'china'
    source = read_image(Path('./Input/' + imagename + '.png'))
    
    write_image('source_' + imagename + '1.png', source)
    
    # Start painting impressionist
    Rs = [8,4,2]
    f_sigma = 0.5
    T = 100 / 255
    f_grid = 1
    f_c= 1
    minLength= 4
    maxLength = 16
    
    canvas = paint(source, Rs, f_sigma, f_grid, T, f_c, minLength, maxLength)
    # blurred = sp.ndimage.gaussian_filter(canvas, [2, 2, 0])
    write_image('Impressionist_' + imagename +'1.png', canvas)
    
    # Start painting Expressionist
    Rs = [8,4,2]
    f_sigma = 0.5
    T = 50 / 255
    f_grid = 1
    f_c= 0.25
    minLength= 10
    maxLength = 16
    Alpha = 0.7
    Jv = 0.3
    
    canvas = paint(source, Rs, f_sigma, f_grid, T, f_c, minLength, maxLength)
    # blurred = sp.ndimage.gaussian_filter(canvas, [2, 2, 0])
    write_image('Expressionist_' + imagename + '1.png', canvas, alpha = Alpha, jv = Jv)
    
    # Start painting Colorist Wash
    Rs = [8,4,2]
    f_sigma = 0.5
    T = 200 / 255
    f_grid = 1
    f_c= 1
    minLength= 4
    maxLength = 16
    Alpha = 0.5
    Jr = 0.3
    Jg = 0.3
    Jb = 0.3
    
    canvas = paint(source, Rs, f_sigma, f_grid, T, f_c, minLength, maxLength)
    # blurred = sp.ndimage.gaussian_filter(canvas, [2, 2, 0])
    write_image('Colorist_' + imagename + '1.png', canvas, alpha = Alpha, jr = Jr, jg = Jg, jb = Jb)
    
    # Start painting Pointillist
    Rs = [4,2]
    f_sigma = 0.5
    T = 100 / 255
    f_grid = 0.5
    f_c= 1
    Alpha = 1.0
    Jv = 1.
    Jh = 0.3
    minLength = 0
    maxLength = 1
    
    canvas = paint(source, Rs, f_sigma, f_grid, T, f_c, minLength, maxLength)
    # blurred = sp.ndimage.gaussian_filter(canvas, [2, 2, 0])
    write_image('Pointillist_' + imagename + '1.png', canvas, alpha = Alpha, jv = Jv, jh = Jh)
    
    
    print('I\'m answering question 1 now!')
    
    ...


def answer_q2() -> None:
    # Load the image
    imagename = 'photo'
    source = read_image(Path('./Input/' + imagename + '.png'))
    
    write_image('source_' + imagename + '.png', source)
    
    # Start painting impressionist
    Rs = [8,4,2]
    f_sigma = 0.5
    T = 20 / 255
    f_grid = 1
    f_c= 1
    minLength= 4
    maxLength = 16
    
    canvas = paint(source, Rs, f_sigma, f_grid, T, f_c, minLength, maxLength)
    # blurred = sp.ndimage.gaussian_filter(canvas, [2, 2, 0])
    write_image('Impressionist_' + imagename +'.png', canvas)
    
    
    # Start painting Pointillist
    Rs = [4,2]
    f_sigma = 0.5
    T = 20 / 255
    f_grid = 0.5
    f_c= 1
    Alpha = 1.0
    Jv = 0.0
    Jh = 0.0
    minLength = 0
    maxLength = 0
    
    canvas = paint(source, Rs, f_sigma, f_grid, T, f_c, minLength, maxLength)
    # blurred = sp.ndimage.gaussian_filter(canvas, [2, 2, 0])
    write_image('Pointillist_' + imagename + '.png', canvas, alpha = Alpha, jv = Jv, jh = Jh)
    print('I\'m answering question 2 now!')
#------------------------------------------------------------------------------
# Define helper functions
def paint(sourceimage, Rs, f_sigma, f_grid, T, f_c, minLength, maxLength):
    canvas = -100 * np.ones(sourceimage.shape)
    print('Canvas shape: ', canvas.shape)
    
    # Use each painter size to paint each layer
    for R in Rs:
        # Apply Gaussian blur
        referenceImage = sp.ndimage.gaussian_filter(sourceimage, [f_sigma * R, f_sigma * R, 0])
        
        # Paint a layer
        paintlayer_stroke(canvas, referenceImage, R, f_grid, T, f_c, minLength, maxLength)
        
    return canvas

# Paint a layer
def paintlayer(canvas, referenceImage, R, f_grid, T):
    
    # Size of the image
    height = canvas.shape[0]
    width = canvas.shape[1]
    
    # List of stroke positions, initially empty
    S = []
    
    # Difference between reference and canvas
    D = np.linalg.norm(canvas - referenceImage, axis = 2)
    
    # Grid window size
    grid = round(f_grid * R)
    
    i = 0
    
    # Find all positions to put on dots
    while i <= height - grid:
        j = 0
        while j <= width - grid:
            areaError = np.sum(D[i:i+grid, j:j+grid]) / grid**2
            if areaError > T:
                
            
                ind = np.unravel_index(np.argmax(D[i:i+grid, j:j+grid], axis=None), D[i:i+grid, j:j+grid].shape)
                ii = ind[0] + i
                jj = ind[1] + j
                S.append([ii,jj])
                # print('Error, T, i, j: ', areaError, T ,i, j)
            j += grid
        i += grid
        #print('i = ', i)
    
      # Make dots in a random order
    random.shuffle(S)
    # print('Number of dots in this layer: ', len(S))
    makeDots(canvas, referenceImage, R, S)
    print('Layer painted: R = ', R)


# Using strokes to paint a layer
def paintlayer_stroke(canvas, referenceImage, R, f_grid, T, f_c, minLength, maxLength):
    
    # Size of the image
    height = canvas.shape[0]
    width = canvas.shape[1]
    
    # List of stroke positions, initially empty
    S = []
    
    # Difference between reference and canvas
    D = np.linalg.norm(canvas - referenceImage, axis = 2)
    
    # Grid window size
    grid = round(f_grid * R)
    
    i = 0
    
    # Find all positions to put on dots
    while i <= height - grid:
        j = 0
        while j <= width - grid:
            areaError = np.sum(D[i:i+grid, j:j+grid]) / grid**2
            if areaError > T:
                
            
                ind = np.unravel_index(np.argmax(D[i:i+grid, j:j+grid], axis=None), D[i:i+grid, j:j+grid].shape)
                ii = ind[0] + i
                jj = ind[1] + j
                S.append([ii,jj])
                # print('Error, T, i, j: ', areaError, T ,i, j)
            j += grid
        i += grid
        #print('i = ', i)
    
      # Make dots in a random order
    random.shuffle(S)
    # print('Number of dots in this layer: ', len(S))
    makeStrokes(canvas, referenceImage, R, S, f_c, minLength, maxLength)
    print('Layer painted: R = ', R)


# Paint strokes with radius R on canvas, using S and reference image
def makeStrokes(canvas, referenceImage, R, S, f_c, minLength, maxLength):
    
    # Size of canvas
    height = canvas.shape[0]
    width = canvas.shape[1]
    
    # Create the masks, decay linearly to 0.7 of initial radius as stroke length increases
    R = int(R)
    
    # List of masks
    masks = []
    
    for k in range(maxLength):
        mask = np.ones([2*R+1, 2*R+1])
        r2 = ((1 - 0.3 * k / maxLength) * R) ** 2
        for i in range(2*R+1):
            for j in range(2*R+1):
                if (i - R)**2 + (j - R)**2 > r2:
                    mask[i,j] = 0
    
        # Stack to have 3 dimensions
        mask = np.dstack((mask, mask, mask))
        masks.append(mask)
    
    
    # Gray of reference for gradient computing
    refGray = 0.3 * referenceImage[:,:,0] + 0.59 * referenceImage[:,:,1] + 0.11 * referenceImage[:,:,2]
    
    grad_ref = np.gradient(refGray)
    
    # Gradient of the layer
    grad_ref = np.dstack((grad_ref[0], grad_ref[1]))
    
    # Paint strokes with start points in s
    for s in S:
    
        # Grab the stroke color
        color = referenceImage[s[0], s[1], :]
        this_stroke = np.dstack((color[0]*np.ones([2 * R + 1, 2 * R + 1]),
        color[1]*np.ones([2 * R + 1, 2 * R + 1]),
        color[2]*np.ones([2 * R + 1, 2 * R + 1])))
        
        # Start point
        x = s[0]
        y = s[1]
        lastDx = 0.
        lastDy = 0.
        points = []
        points.append((x, y))
        
        
        for i in range(maxLength):
            
            # Return comdition
            if (i > minLength) and (np.linalg.norm(referenceImage[x, y, :] - canvas[x,y,:]) < np.linalg.norm(referenceImage[x, y, :] - color)):
                break
            
            # Normalize gradient
            temp = np.linalg.norm(grad_ref[x, y, :])
            if (temp == 0):
                break
            
            # Normal direction
            dy = grad_ref[x,y,0] / temp
            dx = - grad_ref[x,y,1] / temp
            
            if (lastDx * dx + lastDy * dy < 0):
                dx = -dx
                dy = -dy
            
            # Filter stroke direction
            dx = f_c * dx + (1 - f_c) * lastDx
            dy = f_c * dy + (1 - f_c) * lastDy
            temp = np.sqrt(dx**2 + dy**2)
            dx = dx / temp
            dy = dy / temp
            
            x = round(x + R*dx)
            y = round(y + R*dy)
            
            if (x < 0 or x >= height or y < 0 or y >= width):
                break
            
            # Push this point into the list
            points.append((x,y))
            lastDx = dx
            lastDy = dy
        
        # Make this stroke on canvas
        makeAstroke(canvas, R, points, masks, this_stroke)
    
    
# Make a stroke
def makeAstroke(canvas, R, points, masks, this_stroke):
    
    R = int(R)
    for i, p in enumerate(points):
        try:
            canvas[p[0] - R : p[0] + R + 1, p[1] - R : p[1] + R + 1, :] = (1 - masks[i]) * canvas[p[0] - R : p[0] + R + 1, p[1] - R : p[1] + R + 1, :] + masks[i] * this_stroke
            '''
            print('This dot shape: ', this_dot.shape)
            print('Mask shape: ', mask.shape)
            print('Canvas patch shape: ', canvas[s[0] - R : s[0] + R + 1, s[1] - R : s[1] + R + 1, :].shape)
            print('s = ', s)
            '''
        except:
            pass
            

# Paint dots with radius R on canvas using referenceImage and S
def makeDots(canvas, referenceImage, R, S):
    R = int(R)
    mask = np.ones([2*R+1, 2*R+1])
    for i in range(2*R+1):
        for j in range(2*R+1):
            if (i - R)**2 + (j - R)**2 > R**2:
                mask[i,j] = 0
    
    # Stack to have 3 dimensions
    mask = np.dstack((mask, mask, mask))
    
    # print('S = ', S)
    for s in S:
        # Grab the target color
        color = referenceImage[s[0], s[1], :]
        this_dot = np.dstack((color[0]*np.ones([2 * R + 1, 2 * R + 1]),
        color[1]*np.ones([2 * R + 1, 2 * R + 1]),
        color[2]*np.ones([2 * R + 1, 2 * R + 1])))
        
        try:
            canvas[s[0] - R : s[0] + R + 1, s[1] - R : s[1] + R + 1, :] = (1 - mask) * canvas[s[0] - R : s[0] + R + 1, s[1] - R : s[1] + R + 1, :] + mask * this_dot
            '''
            print('This dot shape: ', this_dot.shape)
            print('Mask shape: ', mask.shape)
            print('Canvas patch shape: ', canvas[s[0] - R : s[0] + R + 1, s[1] - R : s[1] + R + 1, :].shape)
            print('s = ', s)
            '''
        except:
            pass
    
    
    
    
#------------------------------------------------------------------------------
# Answer every question if this script is being invoked as the main module.

if __name__ == '__main__':
    for key, value in list(globals().items()):
        if key.startswith('answer_'):
            value()
