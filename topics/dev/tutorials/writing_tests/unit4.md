
> <solution-title></solution-title>
> 
> Factoring out code into a fixture: we declare a fixture using pytest's built-in `fixture` decorator,
> and then simply pass it as an argument to our tests.
> 
> ```python
> @pytest.fixture
> def bytesize():
>     return ByteSize(1_000_000_000_000_000)
> 
> 
> def test_bytesize_to_unit(bytesize):
>     assert bytesize.to_unit('k') == '1000000000000K'
>     assert bytesize.to_unit('m') == '1000000000M'
>     assert bytesize.to_unit('g') == '1000000G'
>     assert bytesize.to_unit('t') == '1000T'
>     assert bytesize.to_unit('p') == '1P'
>     assert bytesize.to_unit('e') == '0E'
> 
>     assert bytesize.to_unit('ki') == '976562500000KI'
>     assert bytesize.to_unit('mi') == '953674316MI'
>     assert bytesize.to_unit('gi') == '931322GI'
>     assert bytesize.to_unit('ti') == '909TI'
>     assert bytesize.to_unit('pi') == '0PI'
>     assert bytesize.to_unit('ei') == '0EI'
> 
>     assert bytesize.to_unit() == '1000000000000000'  # no unit
>     assert bytesize.to_unit('K') == '1000000000000K'  # uppercase unit
> 
> 
> def test_bytesize_to_unit_as_str(bytesize):
>     assert bytesize.to_unit('k', as_string=False) == 1000000000000  # not to_string
> 
> 
> def test_bytesize_to_unit_invalid(bytesize):
>     with pytest.raises(KeyError):
>         bytesize.to_unit('invalid')
> ```
{: .solution }
