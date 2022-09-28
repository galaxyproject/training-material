
> <solution-title></solution-title>
> 
> Factor out the logic for obtaining both lists:
> 1. Create 2 functions that take a transaction object as an argument and return the lists.
> 2. Call the functions from the `validate_email` function's body, assigning the results to
>    local variables.
> 3. Replace all calls to `trans.app.config.[the list]` with the variables you've just assigned.
> 
> 
> 
> ```python
> def get_email_domain_allowlist_content(trans):
>     return trans.app.config.email_domain_allowlist_content
> 
> 
> def get_email_domain_blocklist_content(trans):
>     return trans.app.config.email_domain_blocklist_content
> 
> 
> def validate_email(trans, email, user=None, check_dup=True, allow_empty=False, validate_domain=False):
>     """
>     Validates the email format, also checks whether the domain is blocklisted in the disposable domains configuration.
>     """
>     if (user and user.email == email) or (email == "" and allow_empty):
>         return ""
>     message = validate_email_str(email)
>     if not message and validate_domain:
>         domain = extract_domain(email)
>         message = validate_domain_resolves(domain)
> 
>     if not message and check_dup and check_for_existing_email(trans, email):
>         message = f"User with email '{email}' already exists."
> 
>     email_domain_allowlist_content = get_email_domain_allowlist_content(trans)
>     email_domain_blocklist_content = get_email_domain_blocklist_content(trans)
> 
>     if not message:
>         # If the allowlist is not empty filter out any domain not in the list and ignore blocklist.
>         if email_domain_allowlist_content is not None:
>             domain = extract_domain(email)
>             if domain not in email_domain_allowlist_content:
>                 message = "Please enter an allowed domain email address for this server."
>         # If the blocklist is not empty filter out the disposable domains.
>         elif email_domain_blocklist_content is not None:
>             domain = extract_domain(email, base_only=True)
>             if domain in email_domain_blocklist_content:
>                 message = "Please enter your permanent email address."
> 
>     return message
> ```
{: .solution }
