
> <solution-title></solution-title>
> 
> Final code listing (excluding the 2 added imports).
> 
> ```python
> @pytest.fixture
> def patch_allowlist(monkeypatch):
>     monkeypatch.setattr(validate_user_input, "get_email_domain_allowlist_content", lambda a: None)
> 
> @pytest.fixture
> def patch_blocklist(monkeypatch):
>     monkeypatch.setattr(validate_user_input, "get_email_domain_blocklist_content", lambda a: None)
> 
> @pytest.fixture
> def patch_check_existing(monkeypatch):
>     monkeypatch.setattr(validate_user_input, "check_for_existing_email", lambda a, b: False)
> 
> 
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
> 
>     def test_ignore_blocklist_if_allowlist_not_empty(self, monkeypatch, patch_check_existing):
>         allowlist = ['my_domain.com']
>         monkeypatch.setattr(validate_user_input, "get_email_domain_allowlist_content", lambda a: allowlist)
> 
>         # but add that domain to blocklist too!
>         blocklist = ['my_domain.com']
>         monkeypatch.setattr(validate_user_input, "get_email_domain_blocklist_content", lambda a: blocklist)
> 
>         assert validate_email(None, "email@my_domain.com") == ""  # we expect blocklist to be ignored
> 
>     def test_blocklist_when_allowlist_is_empty(self, monkeypatch, patch_allowlist, patch_check_existing):
>         blocklist = ['bad_domain.com']
>         monkeypatch.setattr(validation_module, "get_email_domain_blocklist_content", lambda a: blocklist)
>         assert "enter your permanent email" in validate_email(None, "email@bad_domain.com")
> ```
{: .solution }
