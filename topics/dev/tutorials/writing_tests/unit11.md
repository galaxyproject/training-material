
> <solution-title></solution-title>
> 
> First take on testing the allowlist.
> 
> ```python
> def test_validate_email_allowlist_not_empty(monkeypatch):
>     allowlist = ['good_domain.com']
>     monkeypatch.setattr(validate_user_input, "get_email_domain_allowlist_content", lambda a: allowlist)
>     monkeypatch.setattr(validate_user_input, "get_email_domain_blocklist_content", lambda a: None)
>     monkeypatch.setattr(validate_user_input, "check_for_existing_email", lambda a, b: False)
> 
>     assert validate_email(None, "email@good_domain.com") == ""
>     assert "enter an allowed domain" in validate_email(None, "email@bad_domain.com")
> ```
{: .solution }
