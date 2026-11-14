
> <solution-title></solution-title>
> 
> Monkeypatching the factored out function. Note the added import.
> 
> ```python
> from galaxy.security import validate_user_input
>
> # [code omitted]
>
> def test_check_for_existing_email(monkeypatch):
>     monkeypatch.setattr(validate_user_input, "check_for_existing_email", lambda a, b: True)
>     result = validate_email(None, "duplicate_email@example.com")
>     assert "exists" in result
> 
>     monkeypatch.setattr(validate_user_input, "check_for_existing_email", lambda a, b: False)
>     result = validate_email(None, "unique_email@example.com")
>     assert result == ""
> ```
{: .solution }

