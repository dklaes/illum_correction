###############
# @file ldac.py
# @author Douglas Applegate & Thomas Erben & Dominik Klaes
# @date 20/02/2014
#
# @brief Utilities to make accessing LDAC cats easier
###############

# HISTORY INFORMATION:
# ====================
#
# 01.09.2010:
# I included treatment of slices through vectorkeys in
# LDACCat.__getitem__
#
# 09.09.2010:
# I made the module more robust against the non-existence
# of necessary libraries
#
# 28.09.2010:
# Bug fix: The from ... import statement must appear on the
# top of the file.
#
# 20.11.2013:
# significant extensions to the code; first proper implementations
# of LDAC tables and LDAC catalogue objects.
#
# 25.11.2013:
# - In the LDACCat class the header of the original catalogue is
#   preserved when storing the catalogue to a file.
# - I add a method to add HISTORY keywords to the image header of
#   the catalogue
#
# 27.11.2013:
# In the LDACTable I replaces the test for an empty LDAC table. Now
# it is done by testing the data element for 'None'. Before it was
# done by the 'size' method on the HDU. It seems that this size element
# of HDU changed from a method to a simple 'int' in different pyfits
# versions. Hence it is useless for code that should be compatible
# to different versions of pyfits.
#
# 20.02.2014
# Dominik Klaes: I implemented adding tables (appending one table to
# the other one).

"""
Wrapper module to work with LDAC catalogues and tables
"""

from __future__ import with_statement

# standard-library includes: 
import sys

# import non-standard modules:
try:
    import pyfits, numpy
except ImportError:
    sys.stderr.write("This script needs the modules 'pyfits' and 'numpy'!\n")
    sys.stderr.write("see http://www.stsci.edu/resources/software_hardware/pyfits\n")
    sys.stderr.write("and/or http://numpy.scipy.org/\n")
    sys.exit(1)


class LDACCat(object):
    """
    Class to represent an LDAC catalogue
    """

    def __init__(self, cat=None):
        """
        An LDAC catalogue can be instantiated either as an empty catalogue
        or with an existing catalogue on disk.

        >>> a = ldac.LDACCat('mag.cat') # reads the catalogue 'mag.cat' into
                                        # the variable 'a'.
        """

        # The LDACCcat object contains a list of LDAC tables.  We
        # internally also keep the header of the PrimaryHDU. It is
        # reused when the catalogue is saved to a file.

        # for an empty catalogue this list is empty:
        self.ldactables = []
        self.header = None
	self.isprimarydatatable = []

        if cat != None:
            # read tables from a catalogue on disk:
            if type(cat) == type("a"):
                hdulist = pyfits.open(cat)

                for hdu in hdulist:
                    if isinstance(hdu, pyfits.PrimaryHDU) == True:
                        self.header = hdu.header
			if hdu.data.size != 0:
				self.ldactables.append(LDACTable(hdu))
			self.isprimarydatatable.append(1)
                    if isinstance(hdu, pyfits.BinTableHDU) == True:
                        self.ldactables.append(LDACTable(hdu))
                        self.isprimarydatatable.append(0)
                    if isinstance(hdu, pyfits.ImageHDU) == True:
                        self.ldactables.append(LDACTable(hdu))
                        self.isprimarydatatable.append(0)

    def __len__(self):
        """
        return the number of LDAC tables in this catalogue

        >>> b = len(a)  # number of LDAC tables in catalogue 'a'.
        """

        return len(self.ldactables)

    def __getitem__(self, tablename):
        """
        returns the named LDAC table. Returns 'None' if the table does
        not exist.

        Example:
        >>> b = a['OBJECTS'] # returns in 'b' the LDAC table with name
                             # 'OBJECTS' from catalogue 'a'
        """

        result = None
        for table in self.ldactables:
            if table.hdu.name == tablename:
                result = table

        return result

    def __setitem__(self, table):
        """
        adds or replaces an LDAC table in this catalogue
        """

        if isinstance(table, LDACTable):
            # check whether a table with table.name exists already:
            exists = False

            for i in xrange(len(self.ldactables)):
                if self.ldactables[i].hdu.name == table.name:
                    self.ldactables[i] = table
                    exists = True

            if exists == False:
                self.ldactables.append(table)

    def tables(self):
        """
        returns the names of the contained LDAC tables

        >>> c = a.tables()  # gives a list of table names in catalogue 'a'
        """
        
        tablenames = []

        for table in self.ldactables:
            tablenames.append(table.hdu.name)

        return tablenames

    def __iter__(self):
        return self.ldactables.__iter__()

    def __contains__(self, tablename):
        """
        check whether a table with name 'tablename' is present
        in this catalogue
        """

        return tablename in self.tables()

    def has_table(self, tablename):
        """
        check whether a table with name 'tablename' is present
        in this catalogue

        >>> c = a.has_table('OBJECTS') # returns 'True' if a table named
                                       # 'OBJECTS' is in catalogue 'a'
        """

        return self.__contains__(tablename)

    def add_history(self, keyvalue):
        """
        add a history keyword to the header of the catalogue

        >>> a.add_history('Catalogue created on 01/02/2013')
        """

        # create an empty header if necessary
        if self.header == None:
            self.header = pyfits.Header()
            
        # just delegate the work to a pyfits method:    
        self.header.add_history('') # empty line for separation from other
                                    # comment/history lines
        self.header.add_history(keyvalue)

    def saveas(self, file, clobber=False):
        """
        save the LDAC catalogue to a file.

        if clobber=True an existing file is overwritten.

        >>> a.saveas('test.cat') # saves LDAC catalogue 'a' with all its
                                 # tables to file 'test.cat'
        """

	if self.isprimarydatatable[0] == 0:
	        primaryHDU = pyfits.PrimaryHDU(header=self.header)
	        hdulist = pyfits.HDUList([primaryHDU])
        	for table in self.ldactables:
	            hdulist.append(table.hdu)
	else:
		primaryHDU = pyfits.PrimaryHDU(data=self.ldactables[0].hdu.data, header=self.ldactables[0].hdu.header)
	        hdulist = pyfits.HDUList([primaryHDU])
        	for table in self.ldactables:
	            if table != self.ldactables[0]:
			    hdulist.append(table.hdu)

        hdulist.writeto(file, clobber=clobber)
        
                
class LDACTable(object):
    """
    Class to represent an LDAC table
    """

    def __init__(self, hdu=None):
        """
        An LDAC table can be instantiated either as am empty table
        or with a pyfits BinaryTable HDU (existing table).
        """

        if hdu == None:
            self.hdu = pyfits.BinTableHDU()
            self.hdu.data = None

            # We make sure that the table has 'some' proper name:
            self.hdu.name = "DEFAULT"
        else:
            self.hdu = hdu
        
    def __len__(self):
        """
        return the number of table entries (objects)
        """

        # 'self.hdu.data' leads to an exception for an empty catalogue.
        # Hence we check for this first:
        if self.hdu.size() == 0:
            return 0
        else:
            return len(self.hdu.data)

    def __getitem__(self, key):
        """
        returns the contents of an existing LDAC key as numpy array

        Example:
        >>> b = a['Xpos'] # store in 'b' the contents (numpy array)
                          # of key 'Xpos' from table 'a'.
        """

        if type(key) == type(5) or \
                type(key) == type(slice(5)):
            return self.hdu.data[key]

        if type(key) == type("a"):
            # we need to deal with slices through vector keys
            # such as 'MAG_APER(2)'
            startind = key.find("(")
            endind = key.find(")")

            if startind > 0 and endind > 0:
                keyname = key[:startind]
                keyindex = int(key[startind + 1:endind]) - 1

                try:
                   return self.hdu.data.field(keyname)[:,keyindex]
                except AttributeError:
                   raise KeyError(key) 
            else:
                try:
                    return self.hdu.data.field(key)
                except AttributeError:
                    raise KeyError(key)

        raise TypeError

    def __setitem__(self, key, val):
        """
        set values of an LDAC table

        a['Xpos'] = b # sets the key 'Xpos' in the table 'a' to the
                      # values in numpy array 'b'. If the key does
                      # not yet exist it is created.
        """
        # VERY uncomplete implementation for the moment!
        # - we only treat scalars for the moment!
        # - we do not check whether the key types match
        #   when an existing key is overwritten

        # sanity checks: the column name must be a string and
        # the value arrays length must match the table data
        # dimension:
        if type(key) == type("a"):
            # The first condition applies to an empty table:
            if self.hdu.data == None or len(val) == self.hdu.data.size:
                # If necessary add a new column to the table
                if self.__contains__(key) == True:
                    # quite some things might go wrong here
                    # (same data type, etc.)

                    # The following construct of '....(key)[:]' ensures
                    # a 'deep' copy of array element which e need here:
                    self.hdu.data.field(key)[:] = val
                else:
                    # determine format for the new column:
                    colformat=""
                    if numpy.issubdtype(val.dtype, float) == True:
                        colformat="1E"
                    
                    if numpy.issubdtype(val.dtype, int) == True:
                        colformat="1I"
                    
                    # now create the new column and create a 'new' table
                    # with the old plus the new column (I did not find a
                    # way to just append a new column to an existing
                    # table!):
                    newcolumn = pyfits.Column(name=key, format=colformat,
                                              array=val)
                    newtabhdu = pyfits.new_table(self.hdu.columns +
                                                 pyfits.ColDefs([newcolumn]))

                    newtabhdu.name = self.hdu.name

                    self.hdu = newtabhdu

        #raise NotImplementedError

    def __delitem__(self, key):
        raise NotImplementedError

    def keys(self):
        """
        returns the names of the keys contained in this table

        >>> b = a.keys() # store a list of keynames of table 'a' in
                         # 'b'.
        """
        return self.hdu.columns.names

    def __iter__(self):
        return self.hdu.data.__iter__()

    def __contains__(self, item):
        return item in self.keys()

    def has_key(self, key):
        """
        tests whether the table contains a certain key.

        >>> b = a.haskey('Xpos') # returns 'True' if table 'a' contains
                                 # a key with name 'Xpos'.
        """
        return self.__contains__(key)

    def filter(self, mask):
        return LDACTable(pyfits.BinTableHDU(data=self.hdu.data[mask],
                                          header=self.hdu.header))

    def setname(self, name):
        """
        set/change the name of the LDAC table.

        >>> a.name = 'TESTTABLE' # set/change the name of the LDAC table
                                 # in 'a' to 'TESTTABLE'.
        """

        self.hdu.name = name

    def saveas(self, file, clobber=False):
        """
        save the LDAC table as a catalogue. The produced
        catalogue will only consist of this table!

        clobber=True overwrites an existing file with the
        new catalogue

        >>> a.saveas('table.cat') # saves the LDAC table in 'a'
                                  # to file 'table.cat'
        """

        self.hdu.writeto(file, clobber=clobber)

    def __add__(a, b):
        """
        Appends table b to table a and returns a LDAC table.
        
        >>> c = a + b   # appends table b to a and saves it
                        # as a LDAC table again
        """
        
        # First check if both tables have the same number of
        # columns:
        if len(a.keys()) != len(b.keys()):
	  print "Tables have not the same number of columns / keywords!"
	  print "First table has " + str(len(a.keys())) + " colums / keywords."
	  print "Second table has " + str(len(b.keys())) + " colums / keywords."
	  return None

        # Now let's check if all kewords from the first table are also
        # present in the second table and also at the same place!
        for i in range(len(a.keys())):
	  if b.has_key(a.keys()[i]) == False:
	    print "Key " + str(a.keys()[i]) + " is not present in the second table!"
	    return None
        
        arows = a.hdu.data.shape[0]
        brows = b.hdu.data.shape[0]
        nrows = arows + brows
        hdu = pyfits.new_table(a.hdu.columns, nrows=nrows)
        print hdu
	
        for i in a.keys():
          hdu.data.field(i)[arows:] = b.hdu.data.field(i)

        hdu.header = a.hdu.header
        hdu.header.update('NAXIS2', nrows)
        hdu.columns = a.hdu.columns
      
        return LDACTable(hdu)

def openObjects(hdulist, table='OBJECTS'):
    tablehdu = None
    for hdu in hdulist:
        # In a regular LDAC catalogue the primary header
        # does not have an EXTNAME keyword and 'hdu.header['EXTNAME']'
        # leads to a KeyError exception which we just ignore:
        try:
            if table == hdu.header['EXTNAME']:
                tablehdu = hdu
        except KeyError:
            pass

    if tablehdu == None:
        print "Table %s not present in catalogue %s" % (table,
                                                        hdulist.filename())
        print "Creating an empty LDAC table"

    return LDACTable(tablehdu)

def openObjectFile(filename, table='OBJECTS'):
    hdulist = pyfits.open(filename)
    if hdulist is None:
        return None

    return openObjects(hdulist, table)
