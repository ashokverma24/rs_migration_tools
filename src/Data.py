#===========================================================================

""": Data module."""

__version__ = "$Revision: #8 $"

#===========================================================================
# Place all imports after here.
#
import pprint
import copy
#
# Place all imports before here.
#===========================================================================


#===========================================================================
class Data( object ):
   """A class that lets you create abitrary data attributes.
   
   Any keyword args passed in to the constructor become data
   attributes.

   # d = Data( a=1, b='2', c=[1,2,3] )
   # print d.a
   # print d.b
   """

   #-----------------------------------------------------------------------
   def __init__( self, **kw ):
      """: Constructor.
      
      = INPUT VARIABLES
      - kw   Dictionary style keyword arguments.
      """
      self.__dict__ = kw
      
   #-----------------------------------------------------------------------
   def __str__( self ):
      ": Convert the object to a string."
      return pprint.pformat( self.__dict__ )

   #-----------------------------------------------------------------------
   def __repr__( self ):
      ": Convert the object to a string."
      return self.__str__()

   #-----------------------------------------------------------------------
   def __len__( self ):
      """: Returns the number of items contained in the specified instance.
      """
      return len( self.__dict__ )

   #-----------------------------------------------------------------------
   def __copy__( self ):
      """: Return a new duplicate instance of this object.
      """
      return Data( **copy.copy( self.__dict__ ) )

   #-----------------------------------------------------------------------
   def __deepcopy__( self, memo ):
      """: Return a deep-copy duplicate of this object.

      = INPUT VARIABLES
      - memo   A dictionary of deep-copied objects.
      """
      return Data( **copy.deepcopy( self.__dict__, memo ) )

   #-----------------------------------------------------------------------
   def __getitem__( self, key ):
      """: Index into the data structure.

      = INPUT VARIABLES
      - key   The dictionary key of the item to get.
      """
      return self.__dict__.__getitem__( key )

   #-----------------------------------------------------------------------
   def __setitem__( self, key, value ):
      """: Index into the data structure and set a value.

      = INPUT VARIABLES
      - key   The dictionary key of the item to set.
      - value The value to set.
      """
      return self.__dict__.__setitem__( key, value )

   #-----------------------------------------------------------------------
   def __iter__( self ):
      ": Iterate over the data structure."
      return self.__dict__.__iter__()

   #-----------------------------------------------------------------------
   def __eq__( self, rhs ):
      ": Equality comparison."
      if not isinstance( rhs, Data ):
         return False

      if len( self ) != len( rhs ):
         return False
      
      for key in self:
         if key not in rhs:
            return False

         if not ( self[key] == rhs[key] ):
            return False

      return True

   #-----------------------------------------------------------------------
   def __ne__( self, rhs ):
      ": Inequality comparison."
      return not self.__eq__( rhs )
   
   #-----------------------------------------------------------------------
   def toDict( self ):
      """: Convert Data object into standard Python dictionary.

      """
      return self.__dict__

   #-----------------------------------------------------------------------
   @staticmethod
   def fromDict( inDict ):
      """: Construct a Data object from a standard Python dictionary.

      = INPUT VARIABLES
      - inDict   The dictionary from which to construct a Data object.
      """
      newData = Data()
      newData.__dict__ = inDict 
      return newData

   #-----------------------------------------------------------------------
   def get( self, key, default=None ):
      """: Return value for the key if found; otherwise, return
      'default'.  If default is not specified, it defaults to None.
      This method never raises a KeyError.

      = INPUT VARIABLES
      - key       The dictionary key of the item to get.
      - default   The value to return if the key does not exist.
      """
      return self.__dict__.get( key, default )

   #-----------------------------------------------------------------------
   def setdefault( self, key, default=None ):
      """: Get a value or set/return the default.

      Returns the value for the key if found; otherwise, inserts the
      key with a value of default and return default.  The input
      default is None if not set.

      = INPUT VARIABLES
      - key       The dictionary key of the item to set.
      - default   The value to return if the key does not exist.
      """
      return self.__dict__.setdefault( key, default )

   #-----------------------------------------------------------------------
   def keys( self ):
      """: Return a list of  possible keys.

      = RETURN VALUE
      - Returns a list of possible key values.
      """
      return self.__dict__.keys()

   #-----------------------------------------------------------------------
   def copy( self, deep = False ):
      """: Return a new duplicate instance of this object.

      = INPUT VARIABLES
      - deep   If set to true, then this will perform a deep copy.

      = RETURN VALUE
      - A copy of this Data object.
      """
      if deep:
         return self.__deepcopy__( {} )
      else:
         return self.__copy__()

   #-----------------------------------------------------------------------
   def update( self, rhs=None, **kwargs ):
      """: Update the data object from another.

      This will update this Data object with the keys from an input
      Data object or dictionary.  This is identical to dict.update().

      """
      if isinstance( rhs, Data ):
         rhs = rhs.__dict__
         
      self.__dict__.update( rhs, **kwargs )
      
   #-----------------------------------------------------------------------

#===========================================================================
