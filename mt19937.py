#!/usr/bin/env python3
#
# mt19937.py - a quick and dirty implementation of the MT19937 PRNG in Python
#
#    Copyright (C) 2020  Tom Liston - email: tom.liston@bad-wolf-sec.com
#                                   - twitter: @tliston
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see [http://www.gnu.org/licenses/].

import random

# this is simply a python implementation of a standard Mersenne Twister PRNG.
# the parameters used, implement the MT19937 variant of the PRNG, based on the
# Mersenne prime 2^19937−1
# see https://en.wikipedia.org/wiki/Mersenne_Twister for a very good explanation
# of the math behind this...

class mt19937():
    u, d = 11, 0xFFFFFFFF
    s, b = 7, 0x9D2C5680
    t, c = 15, 0xEFC60000
    l = 18
    n = 624

    def my_int32(self, x):
        return(x & 0xFFFFFFFF)

    def __init__(self, seed):
        w = 32
        r = 31
        f = 1812433253
        self.m = 397
        self.a = 0x9908B0DF
        self.MT = [0] * self.n
        self.index = self.n + 1
        self.lower_mask = (1 << r) - 1
        self.upper_mask = self.my_int32(~self.lower_mask)
        self.MT[0] = self.my_int32(seed)
        for i in range(1, self.n):
            self.MT[i] = self.my_int32((f * (self.MT[i - 1] ^ (self.MT[i - 1] >> (w - 2))) + i))

    def extract_number(self):
        if self.index >= self.n:
            self.twist()
            self.index = 0
        y = self.MT[self.index]
        # this implements the so-called "tempering matrix"
        # this, functionally, should alter the output to
        # provide a better, higher-dimensional distribution
        # of the most significant bits in the numbers extracted
        y = y ^ ((y >> self.u) & self.d)
        y = y ^ ((y << self.s) & self.b)
        y = y ^ ((y << self.t) & self.c)
        y = y ^ (y >> self.l)
        self.index += 1
        return self.my_int32(y)

    def twist(self):
        for i in range(0, self.n):
            x = (self.MT[i] & self.upper_mask) + (self.MT[(i + 1) % self.n] & self.lower_mask)
            xA = x >> 1
            if(x % 2) != 0:
                xA = xA ^ self.a
            self.MT[i] = self.MT[(i + self.m) % self.n] ^ xA


# so... guess what! while it isn't necessarily obvious, the
# functioning of the tempering matrix are mathematically
# reversible. this function impliments that...
#
# by using this, we can take the output of the MT PRNG, and turn
# it back into the actual values held within the MT[] array itself
# and therefore, we can "clone" the state of the PRNG from "n"
# generated random numbers...
#
# initially, figuring out the math to do this made my brain hurt.
# simplifying it caused even more pain.
# please don't ask me to explain it...
def untemper(y):
    y ^= y >> mt19937.l
    y ^= y << mt19937.t & mt19937.c
    for i in range(7):
        y ^= y << mt19937.s & mt19937.b
    for i in range(3):
        y ^= y >> mt19937.u & mt19937.d
    return y


if __name__ == "__main__":
    # create our own version of an MT19937 PRNG.
    myprng = mt19937(0)
    # fire up Python's built-in PRNG and seed it with the time...
    print("Seeding Python's built-in PRNG with the time...")
    random.seed()
    # generate some random numbers so we can create a random number of random numbers using Python's built-in PRNG
    # so random...
    count1 = random.randint(2000, 10000)
    count2 = random.randint(2000, 10000)
    print("Generating a random number (%i) of random numbers using Python's built-in PRNG..." % (count1))
    print("We do this just to show that this method doesn't depend on being at a particular starting point.")
    for i in range(count1):
        f = random.randrange(0xFFFFFFFF)
    # clone that sucker...
    print("Generating %i random numbers.\nWe'll use those values to create a clone of the current state of Python's built-in PRNG..." % (mt19937.n))
    for i in range(mt19937.n):
        myprng.MT[i] = untemper(random.randrange(0xFFFFFFFF))
    # check to make sure our cloning worked...
    print("Generating a random number (%i) of additional random numbers using Python's built-in PRNG..." % (count2))
    print("Generating those %i random numbers with our clone as well..." % (count2))
    # generate numbers and throw 'em away...
    for i in range(count2):
        f = random.randrange(0xFFFFFFFF)
        f2 = myprng.extract_number()
    print("Now, we'll test the clone...")
    print("\nPython       Our clone")
    for i in range(20):
        r1 = random.randrange(0xFFFFFFFF)
        r2 = myprng.extract_number()
        print("%10.10i - %10.10i (%r)" % (r1, r2, (r1 == r2)))
        assert(r1 == r2)
