"""
This module contains the classes Section and Option to simplify parsing and
validating ini-files.

:Author: Samuel Scherrer
"""

import sys
from copy import deepcopy

__all__ = ["Section", "Option", "MissingSectionError", "MissingOptionError", "ShapeError", "assert_exists"]


class Section():

    def __init__(self, name, options):
        """
        Parameters
        ----------
        name : str
            Name of the section
        options : list
            List of Options
        """
        self.name = name
        self.options = deepcopy(options)
        self.optionnames = [o.name for o in options]

    def parse(self, section_dict):
        """
        Parses and validates options of the section given a dictionary
        containing key-value pairs of the options. If options are not present,
        default values might be set.

        Parameters
        ----------
        section_dict : dictionary
            Dictionary of section values. This can for example be obtained
            using ``configparser``::

                p = configparser.ConfigParser()
                p.read(filename)
                section_dict = dict(p[section_name])

        Returns
        -------
        section_dict : dictionary
            The same dictionary with parsed and validated values
        """
        self.dict = section_dict
        for option in section_dict:
            if option not in self.optionnames:
                print("Warning: Unknown option: {:s} in section {:s}".format(
                    option, self.name), file=sys.stderr
                )
        for option, name in zip(self.options, self.optionnames):
            self.dict[name] = option.parse(self)
        return self.dict

    def __repr__(self):
        repr =  "Section(name={:s},options=".format(self.name)
        for name in self.optionnames:
            repr += name + ','
        repr += ')'
        return repr


class Option():
    """
    Parameters
    ----------
    name : string
        name of the option
    dtype : type
        type of the value of the option, e.g. float, str, int, or list.
        If dtype is list, every entry of the list must have the same type.
    optional : bool, optional (default: False)
        whether the option is optional or not
    default : optional, (default: None)
        if optional=True, this default value will be used if this option is not
        given.
    shape : int or string or list/tuple of ints and/or strings, optional (only
        required for lists, default: None)
        If dtype is ``list``, this is the shape of the list.
        There are several possibilities how to use this option:

        * if ``shape=n`` where n is a nonnegative integer, the value must be
          a list with length ``n``.
        * if ``shape=-1`` the value can have an arbitrary shape.
        * if ``shape="option1"``, the value of this option must have the same
          shape as "option1". This is especially useful if the shape of
          "option1" is set to -1.
        * if ``shape=[2, 3]``, the value must be a list of lists, where the
          outermost list has length 2 and the inner lists all have length 3.
          This also works for more than 2 dimensions.
        * if ``shape=[2, -1, 3]``, the value must be a list of lists of lists.
          The outermost list must again have length 2, the innermost lists must
          have length 3, and the lists at the intermediate level can have any
          length (even different lengths).
        * if ``shape=[2, "option1", 3]``, the values must again be a list of
          lists similar to above, but now the lists at the intermediate level
          must have the same length as "option1".
        * if ``shape=[2, [1, 2], [[3], [3, 3]]]``, the value must be a list
          of lists of lists. The outermost list must again have length 2. The
          two lists it contains have length 1 and length 2. The innermost lists
          all must have length 3.

        It's also possible to only give the innermost value(s), e.g. for a list
        with ``shape=[2, 3]`` only the value ``18``. This will then be expanded
        to ``[[18, 18, 18], [18, 18, 18]]``. Similarly, ``[18, 19, 20]`` would
        be expanded to ``[[18, 19, 20], [18, 19, 20]]``.
        This expansion obviously doesn't work if ``shape=-1``.
        If ``shape=[2, -1, 3]``, expansion is possible if the given value is
        e.g. ``[[1, 2, 3], [1, 2, 3], [1, 2, 3]]``, but not if only ``[1, 2, 3]``
        is given, since the length of the second dimension must be determined
        from the given value.
    datatype: int, float or string, optional (only required for lists, default: None)
        Type of the innermost elements of the option in case the option is a list.
    validation_func: function of one argument that returns a boolean, optional.
        Returns true if the value (for lists or lists of lists this applies to
        all values) is valid.
    alternatives: list of strings
        List of alternative names for this option. If the option is not found
        in the given dictionary while parsing, these alternatives are used if
        they exist. If an alternative is found (searching happens in the
        supplied order), the option is stored under it's standard name.
    """

    def __init__(self, name, dtype, optional=False, default=None, shape=-1,
                 datatype=None, validation_func=lambda x: True, alternatives=[]):
        self.name = name
        self.dtype = dtype
        self.optional = optional
        self.default = default
        self.shape = shape
        self.datatype = datatype
        self.validation_func = validation_func
        if not isinstance(alternatives, list):
            alternatives = [alternatives]
        self.alternatives = alternatives

        # make sure we have shape and datatype for lists
        if self.dtype == list:
            if self.shape is None:
                raise ValueError('Option ' + self.name + ' has type list, but no shape was given.')
            if self.datatype is None:
                raise ValueError('Option ' + self.name + ' has type list, but no datatype for its elements was given.')

    def __repr__(self):
        return "Option(name={:s}, dtype={:s}, optional={:s}, default={:s}, shape={:s}, datatype={:})".format(
            self.name, str(self.dtype), str(self.optional), str(self.default), str(self.shape), str(self.datatype)
        )

    def parse(self, section):
        """
        Parses the option based on it's attributes.
        """
        # try to find alternatives if they exist
        alternatives = deepcopy(self.alternatives)
        while len(alternatives) != 0 and self.name not in section.dict:
            other_name = alternatives.pop(0)
            if other_name in section.dict:
                section.dict[self.name] = section.dict[other_name]
                del section.dict[other_name]
                break
        if not self.optional:
            assert_exists(self.name, section.dict, section.name)
        if self.name not in section.dict:
            return self.default
        else:
            if self.dtype != list:
                if self.dtype == bool:
                    # this is necessary since ``bool("False")`` returns ``True``.
                    value = parse_bool(section, self.name)
                else:
                    value = self.dtype(section.dict[self.name])
                if not self.validation_func(value):
                    raise ValueError('Invalid input for option ' + self.name +
                                    ' in section ' + section.name)
                return value
            else:

                value = parse_list(section.dict[self.name], self.datatype)

                # value validation
                if not all_true(self.validation_func, value):
                    raise ValueError('Invalid input for option ' + self.name +
                                    ' in section ' + section.name)

                shape = deepcopy(self.shape)

                # now we need to get the correct shape
                if shape == -1:
                    # we don't care for the shape of this
                    if not isinstance(value, list):
                        value = [value]
                    return value

                if isinstance(shape, str):
                    # in this case we simply use the shape of the option with this name
                    if shape not in section.dict:
                        raise ValueError(self.name + ' in ' + section.name + ' has an invalid ' +\
                                         'shape because the options whose shape it should have ' +\
                                         'does not exist. Check your option definitions!')
                        raise ShapeError(self.name, section.name)
                    shape = get_shape(section.dict[shape])
                if isinstance(shape, int):
                    shape = [shape]
                # shape is now a list, but it might still contain strings
                for i in range(len(shape)):
                    if isinstance(shape[i], str):
                        shape[i] = len(section.dict[shape[i]])



                # shape is now either a 'flat' shape, i.e. something like [2, 3, 2],
                # or an expanded shape, e.g. [2, [3, 3], [[2, 2, 2],[2, 2, 2]]]
                # if it's flat, it might contain dimensions with -1 that cannot be
                # autoexpanded. We first need to determine the shape of this dimension.
                if is_flat(shape):
                    real_shape = get_shape(value)
                    if isinstance(real_shape, (list, tuple)):
                        # if it's just a single number we can expand it
                        # Here I'm trying to find the flat shape of the value that was
                        # given in the configuration file.
                        flat_shape_value = try_flattening_shape(real_shape)
                        # It might happen that we cannot flatten the shape, in this
                        # case there are negative values remaining in flat_shape_value.
                        # If there are, this means that there is a dimension
                        # containing lists of different lengths.
                        # In any case I will try to replace any -1 in ``shape``
                        # with the value in ``flat_shape_value``.
                        shape = get_positive_shape(shape, flat_shape_value)
                        # Now we do a test for equality of the asserted shape and
                        # the shape of the value found in the config file. Keep in
                        # mind that there might be -1 values left.
                        if flat_shape_value != shape[-len(flat_shape_value):]:
                            raise ShapeError(self.name, section.name)
                        # If there are -1's left we must ensure that the "depth" of
                        # the given value, i.e. the number of dimensions, is higher
                        # than the ``number of dimensions after the value preceding
                        # the first -1`` + 1 .
                        if any(map(lambda x: x == -1, shape)):
                            depth = numdim(value)
                            mindepth = len(shape) - shape.index(-1) + 1
                            if depth < mindepth:
                                raise ValueError('Option ' + self.name + ' in section ' +
                                                section.name + ' can not be expanded!')
                        shape = expand_shape(shape)

                # Now we have an expanded shape, so only two tasks remain:
                # * auto-expansion
                # * shape validation
                value = expand_to_shape(shape, value)
                if not compare_shapes(shape, get_shape(value)):
                    raise ShapeError(self.name, section.name)
                return value



def numdim(l):
    """
    Returns number or dimensions of the list, assuming it has the same depth
    everywhere.
    """
    if not isinstance(l, (list, tuple)):
        return 0
    if not isinstance(l[-1], (list, tuple)):
        return 1
    else:
        return 1 + numdim(l[-1])


def _get_shape(l):
    depth = numdim(l)
    if depth == 0:
        return -1
    if depth == 1:
        return len(l)
    else:
        return [_get_shape(s) for s in l]

def get_shape(l):
    """
    Returns the expanded shape of a nested list, e.g.:

    >>> get_shape(42)
    1
    >>> get_shape([1, 2])
    [2]
    >>> get_shape([[1, 2], [3, 4, 5]])
    [2, [2, 3]]
    >>> get_shape([ [[1, 2], [3, 4]] , [[5, 6, 7], [8]] ])
    [2, [2, 2], [[2, 2], [3, 1]]]
    """
    s = _get_shape(l)
    if s == -1:
        return 1
    shape = [s]
    for i in range(numdim(s)):
        s = _get_shape(s)
        shape.append(s)
    return shape[::-1]

def is_flat(shape):
    return all(map(lambda x: numdim(x) == 0, shape))

def try_flattening_shape(shape):
    depth = numdim(shape)
    flat = [shape[0]]
    for i in range(1, depth):
        flat_list = unroll(shape[i])
        if all_true(lambda x: x == flat_list[0], flat_list):
            flat.append(flat_list[0])
        else:
            flat.append(-1)
    return flat

def get_positive_shape(shape_in_option, shape_in_file):
    """
    Returns a flat shape without -1, replacing the first part with the given
    shape in the option if it's only partially given in the file.
    """
    shape = deepcopy(shape_in_option)
    for i in range(len(shape_in_file)):
        if shape[len(shape)-1-i] == -1:
            shape[len(shape)-1-i] = shape_in_file[len(shape_in_file)-1-i]
    return shape




def unroll(l):
    if numdim(l) == 0:
        return [l]
    if numdim(l) == 1:
        return l
    else:
        unrolled = []
        for li in l:
            unrolled += unroll(li)
        return unrolled


def expand_shape(shape):
    """
    Expands a flat shape to an expanded shape

    >>> expand_shape([2, 3])
    [2, [3, 3]]
    >>> expand_shape([2, 3, 4])
    [2, [3, 3], [[4, 4, 4], [4, 4, 4]]]
    """
    expanded = [shape[0]]
    for i in range(1, len(shape)):
        next = [shape[i]] * shape[i-1]
        for j in range(1, i):
            next = [next] * shape[j]
        expanded.append(next)
    return expanded

def expand_to_shape(shape, value):
    depth = numdim(value)
    if depth != len(shape):
        value = _expand_to_shape(shape[-1-depth], value)
    return value

def _expand_to_shape(shape, value):
    if isinstance(shape, int):
        value = [value] * shape
        return value
    else:
        return [_expand_to_shape(s, value) for s in shape]

def compare_shapes(s1, s2):
    d1 = numdim(s1)
    d2 = numdim(s2)
    if d1 != d2:
        return False
    else:
        l1 = unroll(s1)
        l2 = unroll(s2)
        if len(l1) != len(l2):
            return False
        else:
            for i in range(len(l1)):
                if l1[i] != l2[i] and l1[i] != -1 and l2[i] != -1:
                    return False
            return True



def assert_exists(option, dictionary, section_name):
    if option not in dictionary:
        raise MissingOptionError(option, section_name)

def all_true(func, l):
    if not isinstance(l, (list,tuple)):
        return func(l)
    else:
        return all([all_true(func, li) for li in l])


def parse_list(string, dtype):
    """
    Parses string of form ``'[as, df]'`` to list of type ``dtype``, e.g. for
    ``dtype=str`` it returns `['as', 'df']``. This works also for a list of lists.

    Parameters
    ----------
    string : string
        Has form ``'[some, entries]``. All entries in the list must have the same type.
    dtype : type converter
        Converter to the type the entries should have, e.g. ``str`` or ``int``

    Returns
    -------
    list of elements of type ``dtype`` or list of lists of elements of type
    ``dtype``. Returns a single value if only a single value is given.
    """
    # l =  string.replace('[', '').replace(']', '').replace(' ', '').split(',')
    s = string.replace(' ', '') # remove all spaces first
    if s[0] == '[': # it's not only a single item
        s = s[1:-1] # remove [ and ] from start and end only
    else: # it's just a single item
        return dtype(s)
    if s[0] == '[': # it's a list of lists
        splitted = s.split('],')
        for i in range(len(splitted)-1):
            splitted[i] += ']' # splitting removed the closing bracket from all but the last item
        l = list(map(lambda x: parse_list(x, dtype), splitted))
    else:
        splitted = s.split(',')
        l = list(map(dtype, splitted))
    return l

def parse_bool(section, optionname):
    """
    Parses a string option as bool. Possible options are "True"/"False",
    "yes"/"no", "1"/"0".
    """
    string = section.dict[optionname]
    if string.lower() == "true" or string.lower() == "yes":
        return True
    elif string.lower() == "false" or string.lower() == "no":
        return False
    elif string.isdigit():
        return bool(int(string))
    else:
        raise ValueError("Option " + optionname + " in section " + section.name
                         + " is not a valid boolean!")



class ShapeError(Exception):
    def __init__(self, optionname, sectionname):
        message = "Option " + optionname + " in section " + sectionname + " has wrong shape!"
        super().__init__(message)

class MissingOptionError(Exception):
    def __init__(self, optionname, sectionname):
        message = "Option " + optionname + " in section " + sectionname + " is missing!"
        super().__init__(message)

class MissingSectionError(Exception):
    def __init__(self, sectionname):
        message = "Section " + sectionname + " is missing!"
        super().__init__(message)
