
Changes to file ``lib/galaxy/model/__init__.py``:

```diff
index 873edab..0fd1437 100644
--- a/lib/galaxy/model/__init__.py
+++ b/lib/galaxy/model/__init__.py
@@ -6591,6 +6591,11 @@ class UserPreference(RepresentById):
         self.value = value
 
 
+class UserFavoriteExtension(RepresentById):
+    def __init__(self, value=None):
+        self.value = value
+
+
 class UserAction(RepresentById):
     def __init__(self, id=None, create_time=None, user_id=None, session_id=None, action=None, params=None, context=None):
         self.id = id

```
