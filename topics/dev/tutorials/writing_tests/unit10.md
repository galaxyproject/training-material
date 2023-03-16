
> <solution-title></solution-title>
> 
> Monkeypatch the email lists: return `None` for both, then override as needed:
> 
> ```python
> def test_check_for_existing_email(monkeypatch):
>     monkeypatch.setattr(validate_user_input, "get_email_domain_allowlist_content", lambda a: None)
>     monkeypatch.setattr(validate_user_input, "get_email_domain_blocklist_content", lambda a: None)
> 
>     monkeypatch.setattr(validate_user_input, "check_for_existing_email", lambda a, b: True)
>     result = validate_email(None, "duplicate_email@example.com")
>     assert "exists" in result
> 
>     monkeypatch.setattr(validate_user_input, "check_for_existing_email", lambda a, b: False)
>     result = validate_email(None, "unique_email@example.com")
>     assert result == ""
> ```
{: .solution }

