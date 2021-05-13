
Changes to file ``lib/galaxy/model/migrate/versions/0175_add_user_favorite_extensions.py``:

```diff
new file mode 100644
index 0000000..c9f3715
--- /dev/null
+++ b/lib/galaxy/model/migrate/versions/0175_add_user_favorite_extensions.py
@@ -0,0 +1,37 @@
+"""
+Migration script to add the user_favorite_extension table.
+"""
+
+import datetime
+import logging
+
+from sqlalchemy import Column, DateTime, ForeignKey, Integer, MetaData, Table, Text
+
+now = datetime.datetime.utcnow
+log = logging.getLogger(__name__)
+metadata = MetaData()
+
+APIKeys_table = Table("user_favorite_extension", metadata,
+                      Column("id", Integer, primary_key=True),
+                      Column("user_id", Integer, ForeignKey("galaxy_user.id"), index=True),
+                      Column("value", Text, index=True, unique=True))
+
+
+def upgrade(migrate_engine):
+    metadata.bind = migrate_engine
+    print(__doc__)
+    metadata.reflect()
+    try:
+        APIKeys_table.create()
+    except Exception:
+        log.exception("Creating user_favorite_extension table failed.")
+
+
+def downgrade(migrate_engine):
+    metadata.bind = migrate_engine
+    # Load existing tables
+    metadata.reflect()
+    try:
+        APIKeys_table.drop()
+    except Exception:
+        log.exception("Dropping user_favorite_extension table failed.")

```
