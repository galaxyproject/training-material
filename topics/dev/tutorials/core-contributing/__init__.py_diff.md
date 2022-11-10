
> <solution-title>``lib/galaxy/model/__init__.py``</solution-title>
> 
> Possible changes to file ``lib/galaxy/model/__init__.py``:
> 
> ```diff
> index 76004a716e..c5f2ea79a8 100644
> --- a/lib/galaxy/model/__init__.py
> +++ b/lib/galaxy/model/__init__.py
> @@ -593,6 +593,7 @@ class User(Base, Dictifiable, RepresentById):
>              & not_(Role.name == User.email)  # type: ignore[has-type]
>          ),
>      )
> +    favorite_extensions = relationship("UserFavoriteExtension", back_populates="user")
>  
>      preferences: association_proxy  # defined at the end of this module
>  
> @@ -9998,3 +9999,12 @@ def receive_init(target, args, kwargs):
>          if obj:
>              add_object_to_object_session(target, obj)
>              return  # Once is enough.
> +
> +class UserFavoriteExtension(Base):
> +    __tablename__ = "user_favorite_extension"
> +
> +    id = Column(Integer, primary_key=True)
> +    user_id = Column(ForeignKey("galaxy_user.id"))
> +    value = Column(TEXT)
> +    user = relationship("User", back_populates="favorite_extensions")
> ```
{: .solution }
