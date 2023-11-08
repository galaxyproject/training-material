
> <solution-title></solution-title>
> 
> Factoring out the database call.
> 
> Old version:
> 
> ```python
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
>     # Code to be refactored
>     if (
>         not message
>         and check_dup
>         and trans.sa_session.query(trans.app.model.User)
>         .filter(func.lower(trans.app.model.User.table.c.email) == email.lower())
>         .first()
>     ):
>         message = f"User with email '{email}' already exists."
> 
>     if not message:
>         # If the allowlist is not empty filter out any domain not in the list and ignore blocklist.
>         if trans.app.config.email_domain_allowlist_content is not None:
>             domain = extract_domain(email)
>             if domain not in trans.app.config.email_domain_allowlist_content:
>                 message = "Please enter an allowed domain email address for this server."
>         # If the blocklist is not empty filter out the disposable domains.
>         elif trans.app.config.email_domain_blocklist_content is not None:
>             domain = extract_domain(email, base_only=True)
>             if domain in trans.app.config.email_domain_blocklist_content:
>                 message = "Please enter your permanent email address."
> 
>     return message
> ```
> 
> Refactored version:
> 
> ```python
> def check_for_existing_email(trans, email):
>     user_record = (trans.sa_session.query(trans.app.model.User)
>         .filter(func.lower(trans.app.model.User.table.c.email) == email.lower())
>         .first()
>     )
>     return bool(user_record)
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
>     # Refactored code
>     if not message and check_dup and check_for_existing_email(trans, email):
>         message = f"User with email '{email}' already exists."
> 
>     if not message:
>         # If the allowlist is not empty filter out any domain not in the list and ignore blocklist.
>         if trans.app.config.email_domain_allowlist_content is not None:
>             domain = extract_domain(email)
>             if domain not in trans.app.config.email_domain_allowlist_content:
>                 message = "Please enter an allowed domain email address for this server."
>         # If the blocklist is not empty filter out the disposable domains.
>         elif trans.app.config.email_domain_blocklist_content is not None:
>             domain = extract_domain(email, base_only=True)
>             if domain in trans.app.config.email_domain_blocklist_content:
>                 message = "Please enter your permanent email address."
> 
>     return message
> ```
{: .solution }

