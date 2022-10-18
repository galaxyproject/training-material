
> <solution-title></solution-title>
> 
> Example of tests for the `parse_bytesize` function:
> 
> ```python
> from galaxy.util.bytesize import ByteSize, parse_bytesize
> 
> 
> def test_parse_bytesize_int_or_float():
>     assert parse_bytesize(42) == 42
>     assert parse_bytesize(42.5) == 42.5
> 
> 
> def test_parse_bytesize_str_no_suffix():
>     assert parse_bytesize('42') == 42
>     assert parse_bytesize('42.5') == 42.5
> 
> 
> def test_parse_bytesize_str_suffix():
>     assert parse_bytesize('10K') == 10000
>     assert parse_bytesize('10 KI') == 10240
> ```
{: .solution }
