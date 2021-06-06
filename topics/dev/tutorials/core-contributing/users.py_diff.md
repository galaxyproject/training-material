
> ### {% icon solution %} ``lib/galaxy/managers/users.py``
> 
> Possible changes to file ``lib/galaxy/managers/users.py``:
> 
> ```diff
> index 3407a6bd75..f15bb4c5d6 100644
> --- a/lib/galaxy/managers/users.py
> +++ b/lib/galaxy/managers/users.py
> @@ -590,6 +590,23 @@ class UserManager(base.ModelManager, deletable.PurgableManagerMixin):
>                  log.exception('Subscribing to the mailing list has failed.')
>                  return "Subscribing to the mailing list has failed."
>  
> +    def get_favorite_extensions(self, user):
> +        return [fe.value for fe in user.favorite_extensions]
> +
> +    def add_favorite_extension(self, user, extension):
> +        fe = model.UserFavoriteExtension(value=extension)
> +        user.favorite_extensions.append(fe)
> +        self.session().add(user)
> +        self.session().flush()
> +
> +    def delete_favorite_extension(self, user, extension):
> +        fes = [fe for fe in user.favorite_extensions if fe.value == extension]
> +        if len(fes) == 0:
> +            raise exceptions.RequestParameterInvalidException("Attempted to unfavorite extension not marked as a favorite.")
> +        fe = fes[0]
> +        self.session().delete(fe)
> +        self.session().flush()
> +
>      def activate(self, user):
>          user.active = True
>          self.session().add(user)
> ```
{: .solution }
