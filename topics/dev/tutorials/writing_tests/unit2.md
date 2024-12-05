
> <solution-title></solution-title>
> 
> Example of a test verifying a raised exception. Make sure you add the `pytest` import.
> 
> ```python
> import pytest
> 
> # [code omitted]
> 
> def test_parse_bytesize_str_invalid():
>     with pytest.raises(ValueError):
>         parse_bytesize('10 invalid-suffix')
>     with pytest.raises(ValueError):
>         parse_bytesize('invalid-value')
> ```
{: .solution }
