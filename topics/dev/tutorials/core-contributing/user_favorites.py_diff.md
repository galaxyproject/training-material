
> ### {% icon solution %} ``lib/galaxy/webapps/galaxy/api/user_favorites.py``
> 
> Possible changes to file ``lib/galaxy/webapps/galaxy/api/user_favorites.py``:
> 
> ```diff
> new file mode 100644
> index 0000000000..3587a9056e
> --- /dev/null
> +++ b/lib/galaxy/webapps/galaxy/api/user_favorites.py
> @@ -0,0 +1,66 @@
> +"""
> +API operations allowing clients to determine datatype supported by Galaxy.
> +"""
> +import logging
> +from typing import List
> +
> +from fastapi import Path
> +
> +from galaxy.managers.context import ProvidesUserContext
> +from galaxy.managers.users import UserManager
> +from . import (
> +    depends,
> +    DependsOnTrans,
> +    Router,
> +)
> +
> +log = logging.getLogger(__name__)
> +
> +router = Router(tags=['user_favorites'])
> +
> +ExtensionPath: str = Path(
> +    ...,  # Mark this Path parameter as required
> +    title="Extension",
> +    description="Target file extension for target operation."
> +)
> +
> +
> +@router.cbv
> +class FastAPIUserFavorites:
> +    user_manager: UserManager = depends(UserManager)
> +
> +    @router.get(
> +        '/api/users/current/favorites/extensions',
> +        summary="List user favroite data types",
> +        response_description="List of data types",
> +    )
> +    def index(
> +        self,
> +        trans: ProvidesUserContext = DependsOnTrans,
> +    ) -> List[str]:
> +        """Gets the list of user's favorite data types."""
> +        return self.user_manager.get_favorite_extensions(trans.user)
> +
> +    @router.post(
> +        '/api/users/current/favorites/extensions/{extension}',
> +        summary="Mark an extension as the current user's favorite.",
> +        response_description="The extension.",
> +    )
> +    def create(
> +        self,
> +        extension: str = ExtensionPath,
> +        trans: ProvidesUserContext = DependsOnTrans,
> +    ) -> str:
> +        self.user_manager.add_favorite_extension(trans.user, extension)
> +        return extension
> +
> +    @router.delete(
> +        '/api/users/current/favorites/extensions/{extension}',
> +        summary="Unmark an extension as the current user's favorite.",
> +    )
> +    def delete(
> +        self,
> +        extension: str = ExtensionPath,
> +        trans: ProvidesUserContext = DependsOnTrans,
> +    ) -> str:
> +        self.user_manager.delete_favorite_extension(trans.user, extension)
> ```
{: .solution }
