import os
import configparser as cp
from hyvr.input.option_parsing import Option, Section, ShapeError

def test_configfile():
    filename = os.path.join(os.path.dirname(__file__), 'configtest.ini')

    options = [
        Option('option1', list, optional=False, shape=-1, datatype=float),
        Option('string_shape', list, optional=False, shape='option1', datatype=int),
        Option('nested', list, optional=False, shape=[3, 2], datatype=int),
        Option('negative', list, optional=False, shape=[3, -1, 2], datatype=str),
        Option('expand', list, optional=False, shape='option1', datatype=float),
        Option('float', float, optional=False),
        Option('optional', float, optional=True, default=1),
        Option('missing', float, optional=False, alternatives=['existing1', 'existing2']),
        Option('first_minus_1', list, optional=False, shape=[-1, 'option1'], datatype=float),
        Option('second_minus_1', list, optional=False, shape=['option1', -1], datatype=float),
        Option('a_flag', bool, optional=False),
        Option('another_flag', bool, optional=False),
        Option('yet_another_flag', bool, optional=False),
    ]

    section = Section('test', options)
    p = cp.ConfigParser()
    p.read(filename)
    section_dict = section.parse(dict(p['test']))

    assert section_dict['option1'] == [[4., 5.], [6., 7., 9.]]
    assert section_dict['string_shape'] == [[6, 7], [8, 9, 10]]
    assert section_dict['nested'] == [[2, 2], [3, 3], [4, 5]]
    negative = [[["hi", "my"], ["name", "is"], ["Samuel", "Scherrer"]],
                [["hi", "my"], ["name", "is"], ["Samuel", "Scherrer"]],
                [["hi", "my"], ["name", "is"], ["Samuel", "Scherrer"]]]
    assert section_dict['negative'] == negative
    assert section_dict['expand'] == [[0.1, 0.1], [0.1, 0.1, 0.1]]
    assert section_dict['float'] == 2e-4
    assert section_dict['optional'] == 1.0
    assert section_dict['missing'] == 2.0
    assert section_dict['first_minus_1'] == [[2., 3.], [4., 5.], [5., 2.], [4., 5.]]
    assert section_dict['second_minus_1'] == [[2., 3., 4., 5., 6.], [7., 0.]]
    assert section_dict['a_flag']
    assert not section_dict['another_flag']
    assert not section_dict['yet_another_flag']


    options = [
        Option('wrong_shape', list, shape=[2, 3], datatype=float),
    ]

    error = Section('error1', options)
    try:
        section_dict = error.parse(dict(p['error1']))
    except ShapeError as e:
        thrown = True
        assert str(e) == 'Option wrong_shape in section error1 has wrong shape!'
    assert thrown


    options = [
        Option('no_expansion', list, shape=[2, -1, 2], datatype=float),
    ]
    error = Section('error2', options)
    try:
        section_dict = error.parse(dict(p['error2']))
    except ValueError as e:
        thrown = True
        assert str(e) == 'Option no_expansion in section error2 can not be expanded!'
    assert thrown
