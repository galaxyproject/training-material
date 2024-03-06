
> <solution-title></solution-title>
> 
> Testing the blocklist.
> 
> ```python
> def test_blocklist_when_allowlist_is_empty(self, monkeypatch, patch_allowlist, patch_check_existing):
>     blocklist = ['bad_domain.com']
>     monkeypatch.setattr(validation_module, "get_email_domain_blocklist_content", lambda a: blocklist)
>     assert "enter your permanent email" in validate_email(None, "email@bad_domain.com")
> ```
{: .solution }
