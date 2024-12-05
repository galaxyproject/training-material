
> <solution-title></solution-title>
> 
> Example of parametrizing a test. Replace the previous version of `test_bytesize_to_unit` with the
> following. Note that the test function takes 2 additional arguments.
> 
> ```python
> @pytest.mark.parametrize(
>     "unit, expected_value",
>     [
>         ('k', '1000000000000K'),
>         ('m', '1000000000M'),
>         ('g', '1000000G'),
>         ('t', '1000T'),
>         ('p', '1P'),
>         ('e', '0E'),
>         ('ki', '976562500000KI'),
>         ('mi', '953674316MI'),
>         ('gi', '931322GI'),
>         ('ti', '909TI'),
>         ('pi', '0PI'),
>         ('ei', '0EI'),
>         (None, '1000000000000000'),
>         ('K', '1000000000000K'),
>     ]
> )
> def test_bytesize_to_unit(unit, expected_value, bytesize):
>     assert bytesize.to_unit(unit) == expected_value
> ```
{: .solution }

