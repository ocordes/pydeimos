# value_generator

# written by: Oliver Cordes 2017-10-04
# changed by: Oliver Cordes 2017-10-04

import numpy as np

class VTGenerator_value( object ):
    def __init__( self, valtype, min, max, step ):
        if ( valtype == 'INT' ):
            self.vals = np.arange( int(min), int(max), step=step )
        else:
            self.vals = np.arange( float(min), float(max), step=step )

        print self.vals


    def __iter__( self ):
        return self.vals.__iter__()


class VTGenerator( object ):
  def __init__( self ):
      self.values = {}

  def add_value( self, name, type, min, max, step ):
      vals = VTGenerator_value( type, min, max, step )
      self.values[name] = vals


  def __iter__( self ):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = [tuple(self.values[pool]) for pool in self.values]
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)

  def product( self ):
      return self.__iter__()

  def product_dict( self ):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = [tuple(self.values[pool]) for pool in self.values]
    result = [[]]
    names = self.values.keys()
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield dict( zip( names,prod ))



if __name__ == '__main__' :
    print( 'Hallo')

    vtgen = VTGenerator()

    vtgen.add_value( 'x', 'INT', 0, 10, 2 )
    vtgen.add_value( 'y', 'FLOAT', 0, 20, 2.1 )

    for i in vtgen:
        print( i )

    for i in vtgen.product_dict():
        print( i )
