
'''
This module implements a master class, sequence, thta defines behavior of objects
storing numeric sequences.
It provides sevaral subclasses that generate popular sequences such as Fibonacci, Lucas,
or Mersenne numbers, prime numbers, random sequneces of integers, and many more.
'''

## imports

from warnings import warn
from random import sample
from math import floor, sqrt, acos
from math import degrees as math_deg


## master class that provides slots and methods for various sequences

class sequence():

    '''
    A master class for sequences.
    It provides the sequence and index slots, defines the length method,
    enables iterations, as well as defines accessor methods for the first
    and the last element, and methods for sum
    '''

    def __init__(self, sequence):

        ## the constructor method

        if not type(sequence) is list:
            raise TypeError('sequence has to be a list of numbers')

        if not type(sequence[0]) is float:
            if not type(sequence[0]) is int:
                raise TypeError('sequence has to be a list of numbers')
        
        ## definition of the slots

        self.sequence = sequence
        self.n = len(sequence)

    def __str__(self):

        ## appearance of the object

        print_cap = 'Sequence of n = ' + str(self.n) + ' elements'
        print_body = '\nfirst element = ' + str("{:e}".format(self.first()))
        print_tail = '\nlast element = ' + str("{:e}".format(self.last()))
        
        return print_cap + print_body + print_tail

    def __len__(self):

        return len(self.sequence)

    def __iter__(self):

        ## enables iteration

        self.index = 0

        return self

    def __next__(self):

        ## iterates over a sequence of Fibonacci numbers
        ## the iterator is re-usable

        try:

            result = self.sequence[self.index]

        except IndexError:

            self.index = 0
            
            raise StopIteration
        
        self.index += 1

        return result

    def __add__(self, other):

        ## element-wise addition of two sequence objects,
        ## truncated to the shortest one, or addition od a number
        ## to a sequence

        if (not isinstance(other, sequence)) and (not isinstance(other, float)) and (not isinstance(other, int)):
            
            raise TypeError('operators to be sequence objects or numbers')

        seq1 = self.sequence

        if isinstance(other, sequence): 

            min_len = min((len(self), len(other)))
            max_len = max((len(self), len(other)))

            if max_len > min_len:
                
                warn('Lengths differ: truncating')

            seq2 = other.sequence

            res_seq = [seq1[i] + seq2[i] for i in range(0, min_len)]

        else:

            res_seq = [i + other for i in seq1]

        return sequence(res_seq)

    def __sub__(self, other):

        ## element-wise subtraction of two sequence objects, 
        ## truncated to the shortest one, or subtraction of a number

        if (not isinstance(other, sequence)) and (not isinstance(other, float)) and (not isinstance(other, int)):
            
            raise TypeError('operators to be sequence objects or numbers')

        seq1 = self.sequence

        if isinstance(other, sequence): 

            min_len = min((len(self), len(other)))
            max_len = max((len(self), len(other)))

            if max_len > min_len:
            
                warn('Lengths differ: truncating')

            seq1 = self.sequence
            seq2 = other.sequence

            res_seq = [seq1[i] - seq2[i] for i in range(0, min_len)]

        else:

            res_seq = [i - other for i in seq1]

        return sequence(res_seq)

    def __mul__(self, other):

        ## element-wise multiplication of two sequence objects,
        ## truncated to the shortest one, or division by a numeber

        if (not isinstance(other, sequence)) and (not isinstance(other, float)) and (not isinstance(other, int)):
            
            raise TypeError('operators to be sequence objects or numbers')

        seq1 = self.sequence

        if isinstance(other, sequence): 

            min_len = min((len(self), len(other)))
            max_len = max((len(self), len(other)))

            if max_len > min_len:
                
                warn('Lengths differ: truncating')

            seq1 = self.sequence
            seq2 = other.sequence

            res_seq = [seq1[i] * seq2[i] for i in range(0, min_len)]

        else:

            res_seq = [i * other for i in seq1]

        return sequence(res_seq)

    def __truediv__(self, other):

        ## element-wise division of two sequence objects;
        ## truncated to the shortest one

        if (not isinstance(other, sequence)) and (not isinstance(other, float)) and (not isinstance(other, int)):
            
            raise TypeError('operators to be sequence objects or numbers')

        seq1 = self.sequence

        if isinstance(other, sequence): 

            min_len = min((len(self), len(other)))
            max_len = max((len(self), len(other)))

            if max_len > min_len:
                
                warn('Lengths differ: truncating')

            seq1 = self.sequence
            seq2 = other.sequence

            res_seq = [seq1[i] / seq2[i] for i in range(0, min_len)]

        else:

            res_seq = [i / other for i in seq1]

        return sequence(res_seq)

    def __floordiv__(self, other):

        ## element-wise floor division of two sequence objects, 
        ## truncated to the shortest one, or floor division by a number

        if (not isinstance(other, sequence)) and (not isinstance(other, float)) and (not isinstance(other, int)):
            
            raise TypeError('operators to be sequence objects or numbers')

        seq1 = self.sequence

        if isinstance(other, sequence): 

            min_len = min((len(self), len(other)))
            max_len = max((len(self), len(other)))

            if max_len > min_len:
                
                warn('Lengths differ: truncating')

            seq1 = self.sequence
            seq2 = other.sequence

            res_seq = [seq1[i] // seq2[i] for i in range(0, min_len)]

        else:

            res_seq = [i // other for i in seq1]

        return sequence(res_seq)

    def __mod__(self, other):

        ## element-wise modulo of two sequence objects, 
        ## truncated to the shortest one, or modulo of division by a number

        if (not isinstance(other, sequence)) and (not isinstance(other, float)) and (not isinstance(other, int)):
            
            raise TypeError('operators to be sequence objects or numbers')

        seq1 = self.sequence

        if isinstance(other, sequence): 

            min_len = min((len(self), len(other)))
            max_len = max((len(self), len(other)))

            if max_len > min_len:
                
                warn('Lengths differ: truncating')

            seq1 = self.sequence
            seq2 = other.sequence

            res_seq = [seq1[i] % seq2[i] for i in range(0, min_len)]

        else:

            res_seq = [i % other for i in seq1]

        return sequence(res_seq)

    def __pow__(self, other):

        ## element-wise exponent of two sequence objects, 
        ## truncated to the shortest one, or exponent by a number

        if (not isinstance(other, sequence)) and (not isinstance(other, float)) and (not isinstance(other, int)):
            
            raise TypeError('operators to be sequence objects or numbers')

        seq1 = self.sequence

        if isinstance(other, sequence): 

            min_len = min((len(self), len(other)))
            max_len = max((len(self), len(other)))

            if max_len > min_len:
                
                warn('Lengths differ: truncating')

            seq1 = self.sequence
            seq2 = other.sequence

            res_seq = [seq1[i] ** seq2[i] for i in range(0, min_len)]

        else:

            res_seq = [i ** other for i in seq1]

        return sequence(res_seq)

    def min(self):

        ## the smallest element

        return min(self.sequence)

    def max(self):

        ## the largest element

        return max(self.sequence)

    def last(self):

        ## last element of the sequence

        return self.sequence[self.n - 1]

    def first(self):

        ## first element of the sequence

        return self.sequence[0]

    def sum(self):

        ## sum of sequence elements

        return sum(self.sequence)

    def cumsum(self):

        ## cumulative sum of sequence elements,
        ## a sequence object

        cum_seq = [0] * len(self)

        cum_seq[0] = self.sequence[0]

        rec_sum = cum_seq[0]

        for i in range(1, len(self)):

            rec_sum += self.sequence[i]
            cum_seq[i] = rec_sum

        return sequence(cum_seq)

    def diff(self):

        ## a sequence of differences of the n-th and (n-1)th elements
        ## of the sequence.
        ## Returns a sequence object

        seq = self.sequence

        diff_seq = [seq[i] - seq[i - 1] for i in range(1, len(self))]

        return sequence(diff_seq)

    def product(self):

        ## product of the sequence elements

        prod = 1
        
        for i in self:

            prod *= i

        return prod

    def cumproduct(self):

        ## cumulative product of sequence elements.
        ## returns a sequence object

        cum_seq = [0] * len(self)

        cum_seq[0] = self.sequence[0]

        rec_prod = cum_seq[0]

        for i in range(1, len(self)):

            rec_prod *= self.sequence[i]
            cum_seq[i] = rec_prod

        return sequence(cum_seq)

    def ratio(self):

        ## a sequence of ratios of the n-th and (n-1)th elements
        ## of the sequence.
        ## Returns a sequence object

        seq = self.sequence

        ratio_seq = [seq[i] / seq[i - 1] for i in range(1, len(self))]

        return sequence(ratio_seq)

    def leg(self):

        ## a sequence of [s_0, s_1, ..., s_n] numbers the following form:
        ## s_0 = Null and s_i = sqrt(x_i^2 - x_(i - 1)^2), 
        ## where [x_0, x_1, ..., x_n] is the current sequence

        seq = self.sequence

        leg_seq = [sqrt(seq[i]**2 - seq[i - 1]**2) for i in range(1, len(self))]

        return sequence(leg_seq)

    def angle(self, degrees = False):

        ## a sequence of [s_0, s_1, ..., s_n] numbers the following form:
        ## s_0 = Null and s_i = acos(x_(i - 1)/x_i), 
        ## where [x_0, x_1, ..., x_n] is the current sequence

        seq = self.sequence      

        if degrees:

            theta_seq = [math_deg(acos(seq[i - 1]/seq[i])) for i in range(1, len(self))]

        else:

            theta_seq = [acos(seq[i - 1]/seq[i]) for i in range(1, len(self))]

        return sequence(theta_seq)

    def cos(self):

        ## a sequence of [s_0, s_1, ..., s_n] numbers the following form:
        ## s_0 = Null and s_i = x_(i - 1)/x_i), 
        ## where [x_0, x_1, ..., x_n] is the current sequence

        seq = self.sequence

        cos_seq = [seq[i] / seq[i - 1] for i in range(1, len(self))]

        return sequence(cos_seq)
        
                    
## random sequences of integers and monotonous sequences

class randint(sequence):

    '''
    Generates a random sequence of integers from the first to last range
    '''

    def __init__(self, first, last, n):

        ## the sonstructor method

        ## entry control

        if not isinstance(first, int):

            raise TypeError('first has to be an integer')

        if not isinstance(last, int):

            raise TypeError('last has to be an integer')

        if (not isinstance(n, int)) or n < 1:

            raise TypeError('n has to be a positive integer')

        ## generation of the random sequence

        if last > first:

            int_sequence = [i for i in range(first, last + 1)]

        elif last < first:

            int_sequence = [i for i in range(first, last + 1, -1)]

        else:

            raise ValueError('first and last must not be equal')

        if n > len(int_sequence):

            raise ValueError('n too large: it has to be <' + str(len(int_sequence)))
        
        rand_sequence = sample(int_sequence, n)

        ## filling the slots
        
        self.n = n

        if last > first:

            self.sequence = sorted(rand_sequence)

        else:

            self.sequence = sorted(rand_sequence, reverse = True)

class monotone(sequence):

    '''
    generates a monotonous/uniform numeric sequence, i.e. a sequence that
    consists of one number repeated n-times
    '''

    ## constructor

    def __init__(self, first, n):

        ## entry control

        if(not isinstance(first, float)) and (not isinstance(first, int)):

            raise TypeError('first has to be a float or integer')

        if not isinstance(n, int):

            raise TypeError('n has to be an integer')

        if n < 1:

            raise ValueError('n has to be a positive integer')

        ## geberation of the sequnece

        self.sequence = [first for i in range(0, n)]
        self.n = len(self.sequence)

## prime number generator via Eratostenes Sieve and integer factorization

class prime_max(sequence):

    '''
    Generates a sequence of prime numbers from the range 2 to last
    with the Eratostenes Sieve algorithm
    '''

    def __init__(self, last, verbose = False):

        ## entry control

        if not isinstance(last, int):

            raise TypeError('last has to be a integer')

        if last < 2:

            raise ValueError('last has to be larger or equal to 2')

        ## generation of a screening integer set and a set to store the primes

        int_set = {i for i in range(2, last + 1)}
        
        prime_set = set()

        while len(int_set) > 0:

            current_prime = sorted(list(int_set))[0]

            if verbose: print('Checking k = ' + str(current_prime))

            prime_set.add(current_prime)

            for i in int_set:

                if i % current_prime == 0:

                    int_set = int_set - {i}

        ## filling the slots

        self.n = len(prime_set)
        self.sequence = sorted(list(prime_set))

class factint(sequence):

    '''
    Computes a sequence of prime factors of a gven integer (arithmetic factorization)
    '''

    def __init__(self, k):

        ## entry control

        if not isinstance(k, int):

            raise TypeError('k has to be an integer')

        if k < 1:

            raise ValueError('k has to be a positive integer')

        ## factorization

        if k == 1:

            ## special case: k = 1

            self.sequence = [1]

        else:

            ### iteratively genetaing a sequence of primes that may serve as
            ### factors: they are not larger than sqrt(k) + 1

            prime_ceiling = floor(sqrt(k)) + 1

            prime_sequence = prime_max(prime_ceiling)
            prime_list = prime_sequence.sequence
          
            running_k = k
            factor_list = []

            while running_k not in prime_list:

                for i in prime_list:

                    if running_k % i == 0:

                        factor_list.append(i)
                        running_k = running_k / i
                        break

                if i == prime_sequence.last(): break

            factor_list.append(running_k)

            self.sequence = factor_list     


## classes for popular numeric sequences

## Fibonacci numbers

class fibonacci(sequence):

    '''generates a Fibonacci sequence of n length. Starts from 1.'''

    def __init__(self, n):

        ## the constructor method

        ## entry control

        if not type(n) is int:
            raise TypeError('n has to be an integer')

        if n < 1:
            raise ValueError('n has to be larger than 0')

        ## calculating the n-length Fibonacci sequence

        if n in [1, 2]:

            f_sequence = [1, 1]
            f_sequence = f_sequence[0:n-1]            

        else:

            f_sequence = [1] * n
            f_sequence[1] = 1

            for i in range(2, n):

                f_sequence[i] = f_sequence[i - 1] + f_sequence[i - 2]
        
        ## filling the slots
        
        self.n = n
        self.sequence = f_sequence

class lucas(sequence):

    '''Generates a sequence of Lucas numbers. Starts from 1.'''

    def __init__(self, n):

        ## the constructor method

        ## entry control

        if not type(n) is int:
            raise TypeError('n has to be an integer')

        if n < 1:
            raise ValueError('n has to be larger than 0')

        ## calculating the n-length Lucas sequence

        if n in [1, 2]:

            l_sequence = [2, 1]
            l_sequence = l_sequence[0:n - 1]

        else:

            l_sequence = [0] * n
            l_sequence[0] = 2
            l_sequence[1] = 1

            for i in range(2, n):

                l_sequence[i] = l_sequence[i - 1] + l_sequence[i - 2]
        
        ## filling the slots
        
        self.n = n
        self.sequence = l_sequence

class mersenne(sequence):

    '''Generates Mersenne numbers'''

    def __init__(self, n):

        ## the constructor method

        ## entry control

        if not type(n) is int:
            raise TypeError('n has to be an integer')

        if n < 0:
            raise ValueError('n has to be larger or equal than 0')

        ## calculating the n-length Mersenne sequence

        if(n == 0):

            m_sequence = [0]

        else:

            m_sequence = [pow(2, i + 1) - 1 for i in range(0, n)]
        
        ## filling the slots
        
        self.n = n
        self.sequence = m_sequence

class arithmetic(sequence):

    '''
    Generates an arthihetic sequence with a given
    starting number, and increment, and the requested length.
    '''

    def __init__(self, first, inc, n):

        ## constructor method

        ## the constructor method

        ## entry control

        if not type(n) is int:
            raise TypeError('n has to be an integer')

        if n < 1:
            raise ValueError('n has to be larger than 0')

        if(n == 1):

            sequence = [first]

        else:
            
            sequence = [0] * n

            sequence[0] = first

            for i in range(1, n):

                sequence[i] = sequence[i - 1] + inc

        ## filling the slots
        
        self.n = n
        self.sequence = sequence

class geometric(sequence):

    '''
    Generates a geometric sequence with a given
    starting number, a ratio, and the requested length.
    '''

    def __init__(self, first, ratio, n):

        ## the constructor method

        ## entry control

        if not type(n) is int:
            raise TypeError('n has to be an integer')

        if n < 1:
            raise ValueError('n has to be larger than 0')

        if(n == 1):

            sequence = [first]

        else:
            
            sequence = [0] * n

            sequence[0] = first

            for i in range(1, n):

                sequence[i] = sequence[i - 1] * ratio

        ## filling the slots
        
        self.n = n
        self.sequence = sequence

class rootint(sequence):

    '''
    Creates a sequence of roots of a given base for increasing natural integers.
    '''

    def __init__(self, base = 2, n = 10):

        ## the constructor method

        ## entry control

        if not type(n) is int:
            raise TypeError('n has to be an integer')

        if n < 1:
            raise ValueError('n has to be larger than 0')

        if(not isinstance(base, int)) and (not isinstance(base, float)):

            raise TypeError('base has to be float or integer')

        ## generating the sequnece

        root_seq = arithmetic(1, 1, n) ** (1/base)

        self.n = n
        self.sequence = root_seq.sequence

class powint(sequence):

    '''
    Creates a sequence of exponents of a given base for increasing natural integers.
    '''

    def __init__(self, base = 2, n = 10):

        ## the constructor method

        ## entry control

        if not type(n) is int:
            raise TypeError('n has to be an integer')

        if n < 1:
            raise ValueError('n has to be larger than 0')

        if(not isinstance(base, int)) and (not isinstance(base, float)):

            raise TypeError('base has to be float or integer')

        ## generating the sequnece

        pow_seq = arithmetic(1, 1, n) ** base

        self.n = n
        self.sequence = pow_seq.sequence

        

        

        

        
    


    

    

    

    

        

