'''
#ach test needs its own function defn
make use of `assert` to test the results of operations

2 tests per function (e.g. test with expected input and unexpected input)
is an ok way to start to make sure your functions handle unexpected data types
'''
import pytest # only need this for pytest.raises()

def my_fancy_function(y):
    return y ** 2

def test_my_fancy_function_good():
    assert my_fancy_function(5) == 25

def test_my_fancy_function_bad():
    '''this is a stupid test but just to show an example'''
    with pytest.raises(TypeError):
        my_fancy_function('G')
