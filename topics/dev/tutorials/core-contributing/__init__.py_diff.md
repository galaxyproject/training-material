
> ### {% icon solution %} lib/galaxy/model/__init__.py
> 
> Possible changes to file ``lib/galaxy/model/__init__.py``:
> 
> ```diff
> index 2c0b8a4dc4..69a74b0ad2 100644
> --- a/lib/galaxy/model/__init__.py
> +++ b/lib/galaxy/model/__init__.py
> @@ -6617,6 +6617,11 @@ class UserPreference(RepresentById):
>          self.value = value
>  
>  
> +class UserFavoriteExtension(RepresentById):
> +    def __init__(self, value=None):
> +        self.value = value
> +
> +
>  class UserAction(RepresentById):
>      def __init__(self, id=None, create_time=None, user_id=None, session_id=None, action=None, params=None, context=None):
>          self.id = id
> ```
{. :solution }
