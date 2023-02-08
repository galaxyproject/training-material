
> <solution-title></solution-title>
> 
> Group related tests in one class. Note that you need to add `self` as the first
> argument to each test function. You also need to add `self` to the call to `MockUser()`. Leave
> the fixtures outside the class body, so that they can be reused by other test code. 
> 
> ```python
> class TestValidateEmail:
> 
>     class MockUser:
>         pass
>     
>     def test_email_is_empty(self):
>         assert validate_email(None, "", allow_empty=True) == ""
>     
>     def test_email_is_user_email(self):
>         my_email = "foo"
>         my_user = self.MockUser()
>         my_user.email = my_email
>         assert validate_email(None, my_email, user=my_user) == ""
>     
>     def test_check_for_existing_email(self, patch_allowlist, patch_blocklist, monkeypatch):
>         monkeypatch.setattr(validate_user_input, "check_for_existing_email", lambda a, b: True)
>         result = validate_email(None, "duplicate_email@example.com")
>         assert "exists" in result
>     
>         monkeypatch.setattr(validate_user_input, "check_for_existing_email", lambda a, b: False)
>         result = validate_email(None, "unique_email@example.com")
>         assert result == ""
>     
>     def test_validate_email_allowlist_not_empty(self, monkeypatch, patch_check_existing, patch_blocklist):
>         allowlist = ['good_domain.com']
>         monkeypatch.setattr(validate_user_input, "get_email_domain_allowlist_content", lambda a: allowlist)
>     
>         assert validate_email(None, "email@good_domain.com") == ""
>         assert "enter an allowed domain" in validate_email(None, "email@bad_domain.com")
> ```
{: .solution }


