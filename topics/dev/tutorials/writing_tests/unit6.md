
> <solution-title></solution-title>
> 
> Example of replacing an object with a test double:
> 
> ```python
> from galaxy.security.validate_user_input import (
>     extract_domain,
>     validate_domain_resolves,
>     validate_email,  # we've added this import
>     validate_email_str,
>     validate_publicname_str,
> )
> 
> 
> class MockUser:
>     pass
> 
> 
> def test_email_is_empty():
>     assert validate_email(None, "", allow_empty=True) == ""
> 
> def test_email_is_user_email():
>     my_email = "foo"
>     my_user = MockUser()
>     my_user.email = my_email
>     assert validate_email(None, my_email, user=my_user) == ""
> ```
{: .solution }
