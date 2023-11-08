
> <solution-title></solution-title>
> 
> Testing ignoring the blocklist.
> 
> ```python
> def test_ignore_blocklist_if_allowlist_not_empty(self, monkeypatch, patch_check_existing):
>     allowlist = ['my_domain.com']
>     monkeypatch.setattr(validate_user_input, "get_email_domain_allowlist_content", lambda a: allowlist)
> 
>     # but add that domain to blocklist too!
>     blocklist = ['my_domain.com']
>     monkeypatch.setattr(validate_user_input, "get_email_domain_blocklist_content", lambda a: blocklist)
> 
>     assert validate_email(None, "email@my_domain.com") == ""  # we expect blocklist to be ignored
> ```
{: .solution }
