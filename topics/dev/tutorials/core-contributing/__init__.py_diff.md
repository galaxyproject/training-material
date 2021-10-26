
> ### {% icon solution %} ``lib/galaxy/model/__init__.py``
> 
> Possible changes to file ``lib/galaxy/model/__init__.py``:
> 
> ```diff
> index 35e9ac1ba1..36a3ea9cc5 100644
> --- a/lib/galaxy/model/__init__.py
> +++ b/lib/galaxy/model/__init__.py
> @@ -512,6 +512,7 @@ class User(Base, Dictifiable, RepresentById):
>          cascade='all, delete-orphan',
>          collection_class=ordering_list('order_index'))
>      _preferences = relationship('UserPreference', collection_class=attribute_mapped_collection('name'))
> +    favorite_extensions=relationship('UserFavoriteExtension')
>      values = relationship('FormValues',
>          primaryjoin=(lambda: User.form_values_id == FormValues.id))  # type: ignore
>      # Add type hint (will this work w/SA?)
> @@ -8557,6 +8558,12 @@ class UserPreference(Base, RepresentById):
>          self.name = name
>          self.value = value
>  
> +class UserFavoriteExtension(Base, RepresentById):
> +    __tablename__ = 'user_favorite_extension'
> +
> +    id = Column(Integer, primary_key=True)
> +    user_id = Column(Integer, ForeignKey("galaxy_user.id"), index=True)
> +    value = Column(Text)
>  
>  class UserAction(Base, RepresentById):
>      __tablename__ = 'user_action'
> ```
{: .solution }