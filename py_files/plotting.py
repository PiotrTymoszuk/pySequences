'''
Plotting representations of sequences with Turtle Graphics
'''

## imports

from turtle import *
from py_files.classes import *

## segment plots

def segmentPlot(object,
                colors = ('red', 'blue'),
                axis_color = 'gray',
                scales = (1, 1),
                size = None, 
                return_home = True):

    '''
    draws a segment plot: sequence elements are represented
    by vertical lines, the line values correspond to the element
    values, element signs are color coded.
    '''

    if not isinstance(object, sequence):

        raise TypeError('object has to be an instance of sequence class')

    assert isinstance(colors, tuple)
    assert isinstance(axis_color, str)
    assert isinstance(scales, tuple)
    assert isinstance(return_home, bool)

    if size is not None:
        
        if (not isinstance(size, float)) and (not isinstance(size, int)):

            raise TypeError('size has to be float or integer')

    def_pensize = pensize()

    ## plotting the axis

    seq_len = len(object)

    up()
    setheading(0)
    backward((seq_len // 2) * scales[0])

    axis_start = pos()
    
    setheading(0)
    color(axis_color)
    down()
    forward(seq_len * scales[0])
    up()
    setposition(axis_start)

    ## plotting the segments

    for i in object:

        pensize(size)

        if(i == 0):

            color('black')

        elif(i > 0):

            color(colors[0])

        else:

            color(colors[1])

        down()
        setheading(90)

        curr_pos = pos()

        forward(i * scales[1])
        up()
        setposition(curr_pos[0] + 1 * scales[0], curr_pos[1])

    pensize(def_pensize)
    color('black')

    if(return_home):
            up()
            home()
            down()
            
## Theodorus snail

def theodorusSnail(object,
                   direction = 'left',
                   scale = 1.0,
                   size = None,
                   return_home = True):

    '''
    Draws a generalized form of Theodorus Snail for
    a sequence. The genuine snail is obtained for an increasing sequence
    of square roots of following natural integers
    '''

    ## entry control
    
    if not isinstance(object, sequence):

        raise TypeError('object has to be an instance of sequence class')

    if direction not in ('left', 'right'):

        raise ValueError('direction has to be left or right')

    assert isinstance(scale, float) or isinstance(scale, int)
    assert isinstance(return_home, bool)

    if size is not None:
        
        if (not isinstance(size, float)) and (not isinstance(size, int)):

            raise TypeError('size has to be float or integer')

    ## plotting globals

    def_pensize = pensize()
    pensize(size)
    home_position = pos()

    angle_seq = object.angle(degrees = True).sequence
    len_seq = object.sequence
    n_steps = len(angle_seq)

    if direction == 'left':

        turn_fun = left

    else:

        turn_fun = right

    ## plotting

    setheading(90)
    step_position = pos()
    step_angle = 0
    
    for i in range(0, n_steps):

        turn_fun(step_angle)
        down()
        forward(len_seq[i] * scale)
        step_position = pos()
        up()
        setposition(home_position)

        turn_fun(angle_seq[i])
        down()
        forward(len_seq[i+1] * scale)
        setposition(step_position)
        up()

        setposition(home_position)
        setheading(90)
        step_angle += angle_seq[i]
        step_angle = step_angle % 360

    ## tidying up

    pensize(def_pensize)
    color('black')

    if(return_home):
        up()
        home()
        down()   


## spirals

def squareSpirals(object,
                  direction = 'left',
                  scale = 1.0,
                  size = None, 
                  return_home = True):

    '''
    draws a square spiral for a sequence object:
    lengths of the segments correspond to the subsequent elements.
    '''

    if not isinstance(object, sequence):

        raise TypeError('object has to be an instance of sequence class')

    if direction not in ('left', 'right'):

        raise ValueError('direction has to be left or right')

    assert isinstance(scale, float) or isinstance(scale, int)
    assert isinstance(return_home, bool)

    if size is not None:
        
        if (not isinstance(size, float)) and (not isinstance(size, int)):

            raise TypeError('size has to be float or integer')
        

    def_pensize = pensize()
    pensize(size)

    ## plotting

    if direction == 'left':

        turn_fun = left

    else:

        turn_fun = right

    for i in object:

        forward(i * scale)
        turn_fun(90)

    pensize(def_pensize)
    color('black')

    if(return_home):
        up()
        home()
        down()

## Ulam spirals

def ulamSpiral(object,
               direction = 'left',
               scale = 1.0,
               size = None,
               highlight = 'black', 
               return_home = True):

    '''
    draws a square spiral genuinely proposed by Ulam for
    prime numbers. Single steps (1, 2, ...) are invisible,
    values of sequence elements are highlighted by color.
    '''

    ## entry control

    if not isinstance(object, sequence):

        raise TypeError('object has to be an instance of sequence class')

    if direction not in ('left', 'right'):

        raise ValueError('direction has to be left or right')

    assert isinstance(scale, float) or isinstance(scale, int)
    assert isinstance(return_home, bool)
    assert isinstance(highlight, str)

    if size is not None:
        
        if (not isinstance(size, float)) and (not isinstance(size, int)):

            raise TypeError('size has to be float or integer')

    if direction == 'left':

        turn_fun = left

    else:

        turn_fun = right

    def_pensize = pensize()
    pensize(size)

    ## plotting

    segment_len = 0
    counter = 0

    color(highlight)

    while counter < object.last():

        up()
        
        segment_len += 1

        for i in range(0, segment_len):

            if counter in object.sequence: down()

            forward(scale)
            up()
            counter += 1

        turn_fun(90)

    pensize(def_pensize)
    color('black')

    if(return_home):
        up()
        home()
        down()
        

            

            

    

