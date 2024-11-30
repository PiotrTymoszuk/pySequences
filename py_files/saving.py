'''
This script provides tools for saving
turtle graphics as files.
'''

## imports

from turtle import *
from PIL import Image
from PIL import EpsImagePlugin

## functions

def saveImage(path,
              extension = 'png',
              width = 1000,
              height = 1000,
              hide_turtle = True):

    '''
    Saves the current turtle screen as an EPS file
    and a bitmap
    '''

    ## entry control

    assert isinstance(path, str)
    assert isinstance(extension, str)
    assert isinstance(hide_turtle, bool)

    if(not isinstance(width, float)) and (not isinstance(width, int)):

        raise TypeError('width has to be foat or integer')

    if(not isinstance(height, float)) and (not isinstance(height, int)):

        raise TypeError('height has to be foat or integer')

    ## globals and plugins

    if hide_turtle: hideturtle()    

    ## saving an EPS file

    current_screen = getscreen()
    current_canvas = current_screen.getcanvas()
    current_canvas.postscript(file= path + '.eps',
                              width = width,
                              height = height)

    ## saving a bitmap
    
    img = Image.open(path + '.eps') 
    img.save(path + '.' + extension)
    
