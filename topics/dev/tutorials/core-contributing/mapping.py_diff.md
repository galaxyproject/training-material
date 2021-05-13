
Changes to file ``lib/galaxy/model/mapping.py``:

```diff
index e8fc136..9bbfe6d 100644
--- a/lib/galaxy/model/mapping.py
+++ b/lib/galaxy/model/mapping.py
@@ -1573,6 +1573,12 @@ model.UserPreference.table = Table(
     Column("name", Unicode(255), index=True),
     Column("value", Text))
 
+model.UserFavoriteExtension.table = Table(
+    "user_favorite_extension", metadata,
+    Column("id", Integer, primary_key=True),
+    Column("user_id", Integer, ForeignKey("galaxy_user.id"), index=True),
+    Column("value", Text))
+
 model.UserAction.table = Table(
     "user_action", metadata,
     Column("id", Integer, primary_key=True),
@@ -1661,6 +1667,7 @@ simple_mapping(model.WorkerProcess)
 
 # User tables.
 mapper(model.UserPreference, model.UserPreference.table, properties={})
+mapper(model.UserFavoriteExtension, model.UserFavoriteExtension.table, properties={})
 mapper(model.UserAction, model.UserAction.table, properties=dict(
     # user=relation( model.User.mapper )
     user=relation(model.User)
@@ -1927,6 +1934,9 @@ mapper(model.User, model.User.table, properties=dict(
     _preferences=relation(model.UserPreference,
         backref="user",
         collection_class=attribute_mapped_collection('name')),
+    favorite_extensions=relation(model.UserFavoriteExtension,
+        backref="user"
+    ),
     # addresses=relation( UserAddress,
     #     primaryjoin=( User.table.c.id == UserAddress.table.c.user_id ) ),
     values=relation(model.FormValues,

```
