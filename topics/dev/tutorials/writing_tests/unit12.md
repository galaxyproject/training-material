
> <solution-title></solution-title>
> 
> Factor out all common monkeypatching into fixtures. Note the `pytest` import.
> 
> ```python
> import pytest
> 
> # [code omitted]
> 
> @pytest.fixture
> def patch_allowlist(monkeypatch):
>     monkeypatch.setattr(validate_user_input, "get_email_domain_allowlist_content", lambda a: None)
> 
> 
> @pytest.fixture
> def patch_blocklist(monkeypatch):
>     monkeypatch.setattr(validate_user_input, "get_email_domain_blocklist_content", lambda a: None)
> 
> 
> @pytest.fixture
> def patch_check_existing(monkeypatch):
>     monkeypatch.setattr(validate_user_input, "check_for_existing_email", lambda a, b: False)
> 
> 
> def test_email_is_empty():
>     assert validate_email(None, "", allow_empty=True) == ""
> 
> 
> def test_email_is_user_email():
>     my_email = "foo"
>     my_user = MockUser()
>     my_user.email = my_email
>     assert validate_email(None, my_email, user=my_user) == ""
> 
> 
> def test_check_for_existing_email(patch_allowlist, patch_blocklist, monkeypatch):
>     monkeypatch.setattr(validate_user_input, "check_for_existing_email", lambda a, b: True)
>     result = validate_email(None, "duplicate_email@example.com")
>     assert "exists" in result
> 
>     monkeypatch.setattr(validate_user_input, "check_for_existing_email", lambda a, b: False)
>     result = validate_email(None, "unique_email@example.com")
>     assert result == ""
> 
> 
> def test_validate_email_allowlist_not_empty(monkeypatch, patch_check_existing, patch_blocklist):
>     allowlist = ['good_domain.com']
>     monkeypatch.setattr(validate_user_input, "get_email_domain_allowlist_content", lambda a: allowlist)
> 
>     assert validate_email(None, "email@good_domain.com") == ""
>     assert "enter an allowed domain" in validate_email(None, "email@bad_domain.com")
> ```
{: .solution }
